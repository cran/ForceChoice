#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {

constexpr double kProbFloor = 1e-12;
constexpr double kLogProbFloor = -27.631021115928547; // log(1e-12)
constexpr double kMaxLogit = 27.63102111592755;       // log((1 - 1e-12) / 1e-12)

inline double logspace_add2(double a, double b) {
  if (!std::isfinite(a)) return b;
  if (!std::isfinite(b)) return a;
  if (a > b) {
    return a + std::log1p(std::exp(b - a));
  }
  return b + std::log1p(std::exp(a - b));
}

inline double logspace_sum(const std::vector<double>& x) {
  double out = -std::numeric_limits<double>::infinity();
  for (double value : x) {
    out = logspace_add2(out, value);
  }
  return out;
}

inline double inv_logit_stable(double x) {
  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(x);
  return z / (1.0 + z);
}

inline int pattern_item_index(double value, int n_items) {
  if (!std::isfinite(value)) {
    stop("Pattern item index cannot be NA or non-finite.");
  }
  int idx = static_cast<int>(value) - 1;
  if (idx < 0 || idx >= n_items) {
    stop("Pattern item index is outside the range of 'par'.");
  }
  return idx;
}

struct FcBlock {
  int total_rows = 0;
  int total_cols = 0;
  int obs_rows = 0;
  int obs_cols = 0;
  std::vector<int> total;
  std::vector<int> obs;
  std::vector< std::vector<int> > matches;
};

bool full_row_matches(const FcBlock& block, int obs_row, int total_row) {
  if (block.obs_cols == block.total_cols) {
    for (int c = 0; c < block.total_cols; ++c) {
      if (block.obs[obs_row * block.obs_cols + c] !=
          block.total[total_row * block.total_cols + c]) {
        return false;
      }
    }
    return true;
  }
  if (block.obs_cols == 2) {
    return block.obs[obs_row * block.obs_cols] ==
             block.total[total_row * block.total_cols] &&
           block.obs[obs_row * block.obs_cols + 1] ==
             block.total[total_row * block.total_cols + block.total_cols - 1];
  }
  if (block.obs_cols == 1) {
    return block.obs[obs_row * block.obs_cols] ==
           block.total[total_row * block.total_cols];
  }
  stop("Unsupported forced-choice pattern length.");
}

FcBlock parse_block(SEXP patterns_total_sexp, SEXP patterns_sexp, int n_items) {
  NumericMatrix patterns_total = as<NumericMatrix>(patterns_total_sexp);
  NumericMatrix patterns = as<NumericMatrix>(patterns_sexp);

  FcBlock block;
  block.total_rows = patterns_total.nrow();
  block.total_cols = patterns_total.ncol();
  block.obs_rows = patterns.nrow();
  block.obs_cols = patterns.ncol();

  if (block.total_rows <= 0 || block.total_cols <= 0 ||
      block.obs_rows <= 0 || block.obs_cols <= 0) {
    stop("Pattern matrices must have positive dimensions.");
  }
  if (!(block.obs_cols == block.total_cols ||
        block.obs_cols == 2 ||
        block.obs_cols == 1)) {
    stop("Unsupported forced-choice pattern length.");
  }

  block.total.resize(static_cast<std::size_t>(block.total_rows) *
                     static_cast<std::size_t>(block.total_cols));
  block.obs.resize(static_cast<std::size_t>(block.obs_rows) *
                   static_cast<std::size_t>(block.obs_cols));

  for (int r = 0; r < block.total_rows; ++r) {
    for (int c = 0; c < block.total_cols; ++c) {
      block.total[r * block.total_cols + c] =
        pattern_item_index(patterns_total(r, c), n_items);
    }
  }
  for (int r = 0; r < block.obs_rows; ++r) {
    for (int c = 0; c < block.obs_cols; ++c) {
      block.obs[r * block.obs_cols + c] =
        pattern_item_index(patterns(r, c), n_items);
    }
  }

  block.matches.resize(block.obs_rows);
  for (int obs_row = 0; obs_row < block.obs_rows; ++obs_row) {
    for (int total_row = 0; total_row < block.total_rows; ++total_row) {
      if (full_row_matches(block, obs_row, total_row)) {
        block.matches[obs_row].push_back(total_row);
      }
    }
    if (block.obs_cols == block.total_cols && block.matches[obs_row].empty()) {
      stop("'patterns' contains a full ranking absent from 'patterns.total'.");
    }
  }

  return block;
}

inline double fcggum_item_logit(const std::vector<double>& theta,
                                const NumericMatrix& par,
                                int item,
                                int D) {
  double sum_a = 0.0;
  for (int d = 0; d < D; ++d) {
    sum_a += par(item, d);
  }

  const double cum_psi0 = par(item, 2 * D) * sum_a;
  const double cum_psi1 = cum_psi0 + par(item, 2 * D + 1) * sum_a;

  double r_sq = 0.0;
  for (int d = 0; d < D; ++d) {
    const double diff = theta[d] - par(item, D + d);
    const double a_id = par(item, d);
    r_sq += a_id * a_id * diff * diff;
  }
  const double r_i = std::sqrt(r_sq + 1e-10);

  const double log_num0 = logspace_add2(-cum_psi0, 3.0 * r_i - cum_psi0);
  const double log_num1 = logspace_add2(r_i - cum_psi1,
                                        2.0 * r_i - cum_psi1);
  double logit = log_num1 - log_num0;
  if (logit > kMaxLogit) return kMaxLogit;
  if (logit < -kMaxLogit) return -kMaxLogit;
  return logit;
}

void fcggum_logits(const std::vector<double>& theta,
                   const NumericMatrix& par,
                   int D,
                   std::vector<double>& logits) {
  const int I = par.nrow();
  logits.resize(I);
  for (int item = 0; item < I; ++item) {
    logits[item] = fcggum_item_logit(theta, par, item, D);
  }
}

inline double fcmirt_item_logit(const std::vector<double>& theta,
                                const NumericMatrix& par,
                                int item,
                                int D) {
  double eta = -par(item, D);
  for (int d = 0; d < D; ++d) {
    eta += theta[d] * par(item, d);
  }

  const double c = par(item, D + 1);
  const double upper = par(item, D + 2);
  double prob = c + (upper - c) * inv_logit_stable(eta);
  if (prob < kProbFloor) prob = kProbFloor;
  if (prob > 1.0 - kProbFloor) prob = 1.0 - kProbFloor;
  return std::log(prob) - std::log1p(-prob);
}

void fcmirt_logits(const std::vector<double>& theta,
                   const NumericMatrix& par,
                   int D,
                   std::vector<double>& logits) {
  const int I = par.nrow();
  logits.resize(I);
  for (int item = 0; item < I; ++item) {
    logits[item] = fcmirt_item_logit(theta, par, item, D);
  }
}

double block_observed_logprob(const std::vector<double>& logits,
                              const FcBlock& block,
                              int response_index) {
  if (response_index < 0 || response_index >= block.obs_rows) {
    return kLogProbFloor;
  }

  std::vector<double> full_logprob(block.total_rows);
  for (int row = 0; row < block.total_rows; ++row) {
    double row_lp = 0.0;
    for (int pos = 0; pos < block.total_cols - 1; ++pos) {
      const int selected = block.total[row * block.total_cols + pos];

      double max_term = -std::numeric_limits<double>::infinity();
      for (int candidate_pos = pos; candidate_pos < block.total_cols;
           ++candidate_pos) {
        const int candidate =
          block.total[row * block.total_cols + candidate_pos];
        max_term = std::max(max_term, logits[candidate]);
      }

      double denom_sum = 0.0;
      for (int candidate_pos = pos; candidate_pos < block.total_cols;
           ++candidate_pos) {
        const int candidate =
          block.total[row * block.total_cols + candidate_pos];
        denom_sum += std::exp(logits[candidate] - max_term);
      }
      row_lp += logits[selected] - (max_term + std::log(denom_sum));
    }
    full_logprob[row] = row_lp;
  }

  std::vector<double> obs_logprob(block.obs_rows);
  for (int obs_row = 0; obs_row < block.obs_rows; ++obs_row) {
    double lp = -std::numeric_limits<double>::infinity();
    for (int total_row : block.matches[obs_row]) {
      lp = logspace_add2(lp, full_logprob[total_row]);
    }
    obs_logprob[obs_row] = lp;
  }

  const double obs_denom = logspace_sum(obs_logprob);
  double out = obs_logprob[response_index] - obs_denom;
  if (!std::isfinite(out) || out < kLogProbFloor) {
    return kLogProbFloor;
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
double cpp_fcggum_block_loglik_weighted(NumericMatrix theta,
                                        NumericMatrix par_block,
                                        IntegerVector response,
                                        SEXP patterns_total,
                                        SEXP patterns,
                                        NumericVector weight) {
  const int N = theta.nrow();
  const int D = theta.ncol();
  const int I = par_block.nrow();

  if (response.size() != N || weight.size() != N) {
    stop("'response' and 'weight' must have one element per person.");
  }
  if (par_block.ncol() != 2 * D + 2) {
    stop("FCGGUM block parameters must have 2 * ncol(theta) + 2 columns.");
  }

  FcBlock block = parse_block(patterns_total, patterns, I);
  std::vector<double> theta_i(D), logits(I);
  double out = 0.0;

  for (int person = 0; person < N; ++person) {
    for (int d = 0; d < D; ++d) {
      theta_i[d] = theta(person, d);
    }
    fcggum_logits(theta_i, par_block, D, logits);
    const int y = response[person];
    if (y == NA_INTEGER) {
      out += weight[person] * kLogProbFloor;
    } else {
      out += weight[person] * block_observed_logprob(logits, block, y - 1);
    }
  }

  return out;
}

// [[Rcpp::export]]
double cpp_fcmirt_block_loglik_weighted(NumericMatrix theta,
                                        NumericMatrix par_block,
                                        IntegerVector response,
                                        SEXP patterns_total,
                                        SEXP patterns,
                                        NumericVector weight) {
  const int N = theta.nrow();
  const int D = theta.ncol();
  const int I = par_block.nrow();

  if (response.size() != N || weight.size() != N) {
    stop("'response' and 'weight' must have one element per person.");
  }
  if (par_block.ncol() != D + 3) {
    stop("FCMIRT block parameters must have ncol(theta) + 3 columns.");
  }

  FcBlock block = parse_block(patterns_total, patterns, I);
  std::vector<double> theta_i(D), logits(I);
  double out = 0.0;

  for (int person = 0; person < N; ++person) {
    for (int d = 0; d < D; ++d) {
      theta_i[d] = theta(person, d);
    }
    fcmirt_logits(theta_i, par_block, D, logits);
    const int y = response[person];
    if (y == NA_INTEGER) {
      out += weight[person] * kLogProbFloor;
    } else {
      out += weight[person] * block_observed_logprob(logits, block, y - 1);
    }
  }

  return out;
}
