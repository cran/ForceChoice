// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include "istem_gibbs_grid.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {

inline bool is_na_real(double x) {
  return NumericVector::is_na(x);
}

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

inline int item_index(double value, int n_items) {
  int idx = static_cast<int>(value) - 1;
  if (idx < 0 || idx >= n_items) {
    stop("Pattern item index is outside the range of 'par'.");
  }
  return idx;
}

inline double inv_logit_stable(double x) {
  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(x);
  return z / (1.0 + z);
}

std::vector<int> get_length_poly(const NumericMatrix& par, int D) {
  const int I = par.nrow();
  const int max_poly = par.ncol() - 2 * D;
  if (max_poly <= 0) {
    stop("'par' must contain at least one tau column.");
  }

  std::vector<int> length_poly(I);
  for (int i = 0; i < I; ++i) {
    int Ki = 0;
    for (int k = 0; k < max_poly; ++k) {
      if (!is_na_real(par(i, 2 * D + k))) {
        ++Ki;
      }
    }
    if (Ki <= 0) {
      stop("Each item must have at least one non-missing tau value.");
    }
    length_poly[i] = Ki;
  }
  return length_poly;
}

std::vector<int> get_length_poly_mgpcm(const NumericMatrix& par, int D) {
  const int I = par.nrow();
  const int max_poly = par.ncol() - D;
  if (max_poly <= 0) {
    stop("'par' must contain D discrimination columns followed by at least one category intercept column.");
  }

  std::vector<int> length_poly(I);
  for (int i = 0; i < I; ++i) {
    int Ki = 0;
    for (int k = 0; k < max_poly; ++k) {
      if (!is_na_real(par(i, D + k))) {
        ++Ki;
      }
    }
    if (Ki <= 0) {
      stop("Each MGPCM item must have at least one non-missing category intercept.");
    }
    length_poly[i] = Ki;
  }
  return length_poly;
}

NumericMatrix mggum_binary_agree(const NumericMatrix& theta,
                                 const NumericMatrix& par) {
  const int N = theta.nrow();
  const int D = theta.ncol();
  const int I = par.nrow();

  if (par.ncol() != 2 * D + 2) {
    stop("FCGGUM requires binary MGGUM item parameters: ncol(par) must be 2 * ncol(theta) + 2.");
  }

  NumericMatrix agree(N, I);
  for (int i = 0; i < I; ++i) {
    double sum_a = 0.0;
    for (int d = 0; d < D; ++d) {
      sum_a += par(i, d);
    }

    const double cum_psi0 = par(i, 2 * D) * sum_a;
    const double cum_psi1 = cum_psi0 + par(i, 2 * D + 1) * sum_a;

    for (int p = 0; p < N; ++p) {
      double r_sq = 0.0;
      for (int d = 0; d < D; ++d) {
        const double diff = theta(p, d) - par(i, D + d);
        const double a_id = par(i, d);
        r_sq += a_id * a_id * diff * diff;
      }
      const double r_i = std::sqrt(r_sq + 1e-10);

      const double log_num0 = logspace_add2(-cum_psi0, 3.0 * r_i - cum_psi0);
      const double log_num1 = logspace_add2(r_i - cum_psi1, 2.0 * r_i - cum_psi1);
      const double log_denom = logspace_add2(log_num0, log_num1);
      agree(p, i) = std::exp(log_num1 - log_denom);
    }
  }
  return agree;
}

NumericMatrix mirt_binary_agree(const NumericMatrix& theta,
                                const NumericMatrix& par) {
  const int N = theta.nrow();
  const int D = theta.ncol();
  const int I = par.nrow();

  if (D <= 0) {
    stop("'theta' must have at least one column.");
  }
  if (par.ncol() != D + 3) {
    stop("'par' must have D discrimination columns followed by b, c, and d.");
  }

  NumericMatrix agree(N, I);
  for (int p = 0; p < N; ++p) {
    for (int i = 0; i < I; ++i) {
      double eta = -par(i, D);
      for (int d = 0; d < D; ++d) {
        eta += theta(p, d) * par(i, d);
      }
      const double c = par(i, D + 1);
      const double upper = par(i, D + 2);
      agree(p, i) = inv_logit_stable(eta) * (upper - c) + c;
    }
  }
  return agree;
}

NumericMatrix forced_choice_from_agree(const NumericMatrix& agree,
                                       List patterns_total,
                                       List patterns) {
  const int N = agree.nrow();
  const int I_states = agree.ncol();
  const int N_block = patterns_total.size();

  if (patterns.size() != N_block) {
    stop("'patterns.total' and 'patterns' must have the same length.");
  }

  NumericMatrix logit(N, I_states);
  const double eps = 1e-12;
  for (int p = 0; p < N; ++p) {
    for (int i = 0; i < I_states; ++i) {
      double prob_i = agree(p, i);
      if (prob_i < eps) prob_i = eps;
      if (prob_i > 1.0 - eps) prob_i = 1.0 - eps;
      logit(p, i) = std::log(prob_i) - std::log1p(-prob_i);
    }
  }

  int prob_cols = 0;
  for (int b = 0; b < N_block; ++b) {
    NumericMatrix patterns_cur = as<NumericMatrix>(patterns[b]);
    prob_cols += patterns_cur.nrow();
  }

  NumericMatrix prob(N, prob_cols);
  int idx = 0;

  for (int b = 0; b < N_block; ++b) {
    NumericMatrix patterns_total_cur = as<NumericMatrix>(patterns_total[b]);
    NumericMatrix patterns_cur = as<NumericMatrix>(patterns[b]);

    const int bR_total = patterns_total_cur.nrow();
    const int bC_total = patterns_total_cur.ncol();
    const int bR = patterns_cur.nrow();
    const int bC = patterns_cur.ncol();

    if (bR_total <= 0 || bC_total <= 0 || bR <= 0 || bC <= 0) {
      stop("Pattern matrices must have positive dimensions.");
    }

    NumericMatrix log_prob_total(N, bR_total);
    NumericMatrix prob_total(N, bR_total);

    for (int br = 0; br < bR_total; ++br) {
      for (int p = 0; p < N; ++p) {
        double log_prob_br = 0.0;
        for (int bc = 0; bc < bC_total - 1; ++bc) {
          const int selected = item_index(patterns_total_cur(br, bc), I_states);

          double max_term = -std::numeric_limits<double>::infinity();
          for (int candidate_pos = bc; candidate_pos < bC_total; ++candidate_pos) {
            const int candidate = item_index(patterns_total_cur(br, candidate_pos), I_states);
            const double term = logit(p, candidate);
            if (term > max_term) {
              max_term = term;
            }
          }

          double denom_sum = 0.0;
          for (int candidate_pos = bc; candidate_pos < bC_total; ++candidate_pos) {
            const int candidate = item_index(patterns_total_cur(br, candidate_pos), I_states);
            denom_sum += std::exp(logit(p, candidate) - max_term);
          }
          const double log_denom = max_term + std::log(denom_sum);
          log_prob_br += logit(p, selected) - log_denom;
        }
        log_prob_total(p, br) = log_prob_br;
      }
    }

    for (int p = 0; p < N; ++p) {
      double max_log = -std::numeric_limits<double>::infinity();
      for (int br = 0; br < bR_total; ++br) {
        if (log_prob_total(p, br) > max_log) {
          max_log = log_prob_total(p, br);
        }
      }
      double denom_sum = 0.0;
      for (int br = 0; br < bR_total; ++br) {
        denom_sum += std::exp(log_prob_total(p, br) - max_log);
      }
      const double log_denom = max_log + std::log(denom_sum);
      for (int br = 0; br < bR_total; ++br) {
        prob_total(p, br) = std::exp(log_prob_total(p, br) - log_denom);
      }
    }

    if (bC == bC_total) {
      std::vector<int> matches(bR, -1);
      for (int br = 0; br < bR; ++br) {
        for (int fr = 0; fr < bR_total; ++fr) {
          bool same = true;
          for (int c = 0; c < bC_total; ++c) {
            if (item_index(patterns_cur(br, c), I_states) !=
                item_index(patterns_total_cur(fr, c), I_states)) {
              same = false;
              break;
            }
          }
          if (same) {
            matches[br] = fr;
            break;
          }
        }
        if (matches[br] < 0) {
          stop("'patterns' contains a full ranking absent from 'patterns.total'.");
        }
      }
      for (int p = 0; p < N; ++p) {
        for (int br = 0; br < bR; ++br) {
          prob(p, idx + br) = prob_total(p, matches[br]);
        }
      }
    } else if (bC == 2) {
      std::vector< std::vector<int> > matches(bR);
      for (int br = 0; br < bR; ++br) {
        const int first = item_index(patterns_cur(br, 0), I_states);
        const int last = item_index(patterns_cur(br, 1), I_states);
        for (int fr = 0; fr < bR_total; ++fr) {
          const int full_first = item_index(patterns_total_cur(fr, 0), I_states);
          const int full_last = item_index(patterns_total_cur(fr, bC_total - 1), I_states);
          if (full_first == first && full_last == last) {
            matches[br].push_back(fr);
          }
        }
      }
      for (int p = 0; p < N; ++p) {
        for (int br = 0; br < bR; ++br) {
          double s = 0.0;
          for (int fr : matches[br]) {
            s += prob_total(p, fr);
          }
          prob(p, idx + br) = s;
        }
      }
    } else if (bC == 1) {
      std::vector< std::vector<int> > matches(bR);
      for (int br = 0; br < bR; ++br) {
        const int first = item_index(patterns_cur(br, 0), I_states);
        for (int fr = 0; fr < bR_total; ++fr) {
          const int full_first = item_index(patterns_total_cur(fr, 0), I_states);
          if (full_first == first) {
            matches[br].push_back(fr);
          }
        }
      }
      for (int p = 0; p < N; ++p) {
        for (int br = 0; br < bR; ++br) {
          double s = 0.0;
          for (int fr : matches[br]) {
            s += prob_total(p, fr);
          }
          prob(p, idx + br) = s;
        }
      }
    } else {
      stop("Unsupported forced-choice pattern length.");
    }

    for (int p = 0; p < N; ++p) {
      double row_sum = 0.0;
      for (int br = 0; br < bR; ++br) {
        row_sum += prob(p, idx + br);
      }
      if (row_sum <= 0.0 || !std::isfinite(row_sum)) {
        stop("Computed forced-choice probabilities cannot be normalized.");
      }
      for (int br = 0; br < bR; ++br) {
        prob(p, idx + br) /= row_sum;
      }
    }

    idx += bR;
  }

  return prob;
}

constexpr double kFcdcmProbFloor = 1e-12;
constexpr double kFcdcmLogFloor = -27.631021115928547; // log(1e-12)

inline double fcdcm_clip_prob(double x) {
  if (x < kFcdcmProbFloor) return kFcdcmProbFloor;
  if (x > 1.0 - kFcdcmProbFloor) return 1.0 - kFcdcmProbFloor;
  return x;
}

inline int fcdcm_item_index(int value, int n_items) {
  const int idx = value - 1;
  if (idx < 0 || idx >= n_items) {
    stop("FCDCM pattern item index is outside the range of 'zeta'.");
  }
  return idx;
}

inline double fcdcm_block_prob(int zeta_a, int zeta_b,
                               double eta0, double etaAB,
                               double eta_equal0 = 0.5,
                               double eta_equal1 = 0.5) {
  if (zeta_a < zeta_b) return fcdcm_clip_prob(eta0);
  if (zeta_a > zeta_b) return fcdcm_clip_prob(0.5 + etaAB);
  return fcdcm_clip_prob(zeta_a == 0 ? eta_equal0 : eta_equal1);
}

std::vector<int> fcdcm_block_case(const IntegerMatrix& zeta_patterns,
                                  const IntegerMatrix& patterns) {
  const int C = zeta_patterns.nrow();
  const int I = zeta_patterns.ncol();
  const int B = patterns.nrow();
  if (patterns.ncol() != 2) {
    stop("'patterns' must have two columns for FCDCM blocks.");
  }

  std::vector<int> out(static_cast<std::size_t>(C) *
                       static_cast<std::size_t>(B));
  std::vector<int> item_a(B), item_b(B);
  for (int b = 0; b < B; ++b) {
    item_a[b] = fcdcm_item_index(patterns(b, 0), I);
    item_b[b] = fcdcm_item_index(patterns(b, 1), I);
  }

  for (int c = 0; c < C; ++c) {
    for (int b = 0; b < B; ++b) {
      const int za = zeta_patterns(c, item_a[b]);
      const int zb = zeta_patterns(c, item_b[b]);
      out[c * B + b] = (za < zb) ? 0 : ((za > zb) ? 2 : (za == 0 ? 1 : 3));
    }
  }
  return out;
}

void fcdcm_check_dims(const NumericMatrix& theta,
                      const NumericVector& delta1,
                      const NumericVector& delta0,
                      const NumericMatrix& par,
                      const NumericMatrix& alpha_patterns,
                      const IntegerMatrix& zeta_patterns,
                      const IntegerMatrix& patterns) {
  if (theta.ncol() < 1) {
    stop("'theta' must have at least one column.");
  }
  const int D = alpha_patterns.ncol();
  const int C = alpha_patterns.nrow();
  const int B = patterns.nrow();
  if (delta1.size() != D || delta0.size() != D) {
    stop("'delta1' and 'delta0' must have length equal to ncol(alpha_patterns).");
  }
  if (zeta_patterns.nrow() != C) {
    stop("'zeta_patterns' must have the same number of rows as 'alpha_patterns'.");
  }
  if (patterns.ncol() != 2) {
    stop("'patterns' must have two columns for FCDCM blocks.");
  }
  if (par.nrow() != B || (par.ncol() != 2 && par.ncol() != 4)) {
    stop("'par' must be a B x 2 or B x 4 FCDCM parameter matrix.");
  }
}

void fcdcm_class_logprob(const NumericMatrix& par,
                         const std::vector<int>& block_case,
                         int C,
                         int B,
                         std::vector<double>& log_p1,
                         std::vector<double>& log_p0) {
  log_p1.resize(static_cast<std::size_t>(C) * static_cast<std::size_t>(B));
  log_p0.resize(static_cast<std::size_t>(C) * static_cast<std::size_t>(B));
  for (int c = 0; c < C; ++c) {
    for (int b = 0; b < B; ++b) {
      double p = 0.5;
      const int case_id = block_case[c * B + b];
      if (case_id == 0) {
        p = par(b, 0);
      } else if (case_id == 2) {
        p = 0.5 + par(b, 1);
      } else if (par.ncol() == 4 && case_id == 1) {
        p = par(b, 2);
      } else if (par.ncol() == 4 && case_id == 3) {
        p = par(b, 3);
      }
      p = fcdcm_clip_prob(p);
      log_p1[c * B + b] = std::log(p);
      log_p0[c * B + b] = std::log1p(-p);
    }
  }
}

void fcdcm_class_logprob_from_prob(const NumericMatrix& class_prob,
                                   std::vector<double>& log_p1,
                                   std::vector<double>& log_p0) {
  const int C = class_prob.nrow();
  const int B = class_prob.ncol();
  log_p1.resize(static_cast<std::size_t>(C) * static_cast<std::size_t>(B));
  log_p0.resize(static_cast<std::size_t>(C) * static_cast<std::size_t>(B));
  for (int c = 0; c < C; ++c) {
    for (int b = 0; b < B; ++b) {
      const double p = fcdcm_clip_prob(class_prob(c, b));
      log_p1[c * B + b] = std::log(p);
      log_p0[c * B + b] = std::log1p(-p);
    }
  }
}

void fcdcm_class_prob(const NumericMatrix& par,
                      const std::vector<int>& block_case,
                      int C,
                      int B,
                      NumericMatrix& prob) {
  for (int c = 0; c < C; ++c) {
    for (int b = 0; b < B; ++b) {
      double p = 0.5;
      const int case_id = block_case[c * B + b];
      if (case_id == 0) {
        p = par(b, 0);
      } else if (case_id == 2) {
        p = 0.5 + par(b, 1);
      } else if (par.ncol() == 4 && case_id == 1) {
        p = par(b, 2);
      } else if (par.ncol() == 4 && case_id == 3) {
        p = par(b, 3);
      }
      prob(c, b) = fcdcm_clip_prob(p);
    }
  }
}

inline double fcdcm_log_alpha_prob(double theta,
                                   const NumericMatrix& alpha_patterns,
                                   int cls,
                                   const NumericVector& delta1,
                                   const NumericVector& delta0) {
  const int D = alpha_patterns.ncol();
  double lp = 0.0;
  for (int d = 0; d < D; ++d) {
    const double eta = delta1[d] * (theta - delta0[d]);
    if (alpha_patterns(cls, d) > 0.5) {
      lp += std::log(fcdcm_clip_prob(inv_logit_stable(eta)));
    } else {
      lp += std::log(fcdcm_clip_prob(inv_logit_stable(-eta)));
    }
  }
  return lp;
}

double fcdcm_person_loglik(double theta,
                           const IntegerMatrix& response,
                           int person,
                           const NumericVector& delta1,
                           const NumericVector& delta0,
                           const NumericMatrix& alpha_patterns,
                           const std::vector<double>& log_p1,
                           const std::vector<double>& log_p0) {
  const int C = alpha_patterns.nrow();
  const int B = response.ncol();
  std::vector<double> lp(C);
  for (int c = 0; c < C; ++c) {
    double val = fcdcm_log_alpha_prob(theta, alpha_patterns, c,
                                      delta1, delta0);
    for (int b = 0; b < B; ++b) {
      const int y = response(person, b);
      val += (y == 1) ? log_p1[c * B + b] : log_p0[c * B + b];
    }
    lp[c] = val;
  }
  return logspace_sum(lp);
}

void fcdcm_alpha_weights(double theta,
                         const NumericVector& delta1,
                         const NumericVector& delta0,
                         const NumericMatrix& alpha_patterns,
                         std::vector<double>& weight) {
  const int C = alpha_patterns.nrow();
  weight.resize(C);
  double max_lp = -std::numeric_limits<double>::infinity();
  for (int c = 0; c < C; ++c) {
    weight[c] = fcdcm_log_alpha_prob(theta, alpha_patterns, c,
                                     delta1, delta0);
    if (weight[c] > max_lp) max_lp = weight[c];
  }

  double denom = 0.0;
  for (int c = 0; c < C; ++c) {
    weight[c] = std::exp(weight[c] - max_lp);
    denom += weight[c];
  }
  if (denom <= 0.0 || !std::isfinite(denom)) {
    stop("Invalid FCDCM attribute-pattern probabilities.");
  }
  for (int c = 0; c < C; ++c) {
    weight[c] /= denom;
  }
}

double fcdcm_log_alpha_row(double theta,
                           const IntegerMatrix& alpha,
                           int person,
                           const NumericVector& delta1,
                           const NumericVector& delta0) {
  const int D = alpha.ncol();
  double lp = 0.0;
  for (int d = 0; d < D; ++d) {
    const double eta = delta1[d] * (theta - delta0[d]);
    if (alpha(person, d) == 1) {
      lp += std::log(fcdcm_clip_prob(inv_logit_stable(eta)));
    } else {
      lp += std::log(fcdcm_clip_prob(inv_logit_stable(-eta)));
    }
  }
  return lp;
}

} // namespace

// [[Rcpp::export]]
NumericMatrix cpp_model_MIRT(NumericMatrix theta, NumericMatrix par) {
  return mirt_binary_agree(theta, par);
}

// [[Rcpp::export]]
NumericMatrix cpp_model_MGGUM(NumericMatrix theta, NumericMatrix par) {
  const int N = theta.nrow();
  const int D = theta.ncol();
  const int I = par.nrow();

  if (D <= 0) {
    stop("'theta' must have at least one column.");
  }
  if (par.ncol() <= 2 * D) {
    stop("'par' must have D discrimination columns, D delta columns, and tau columns.");
  }

  std::vector<int> length_poly = get_length_poly(par, D);
  int total_cols = 0;
  for (int Ki : length_poly) {
    total_cols += Ki;
  }

  NumericMatrix prob(N, total_cols);
  int idx = 0;
  std::vector<double> cum_psi;
  std::vector<double> log_num;

  for (int i = 0; i < I; ++i) {
    const int Ki = length_poly[i];
    const int Mi = 2 * Ki - 1;

    double sum_a = 0.0;
    for (int d = 0; d < D; ++d) {
      sum_a += par(i, d);
    }

    cum_psi.assign(Ki, 0.0);
    double running_psi = 0.0;
    for (int k = 0; k < Ki; ++k) {
      running_psi += par(i, 2 * D + k) * sum_a;
      cum_psi[k] = running_psi;
    }

    log_num.assign(Ki, 0.0);
    for (int p = 0; p < N; ++p) {
      double r_sq = 0.0;
      for (int d = 0; d < D; ++d) {
        const double diff = theta(p, d) - par(i, D + d);
        const double a_id = par(i, d);
        r_sq += a_id * a_id * diff * diff;
      }
      const double r_i = std::sqrt(r_sq + 1e-10);

      double max_log_num = -std::numeric_limits<double>::infinity();
      for (int k = 0; k < Ki; ++k) {
        const double z = static_cast<double>(k);
        const double arg1 = z * r_i - cum_psi[k];
        const double arg2 = (static_cast<double>(Mi) - z) * r_i - cum_psi[k];
        log_num[k] = logspace_add2(arg1, arg2);
        if (log_num[k] > max_log_num) {
          max_log_num = log_num[k];
        }
      }

      double denom_sum = 0.0;
      for (int k = 0; k < Ki; ++k) {
        denom_sum += std::exp(log_num[k] - max_log_num);
      }
      const double log_denom = max_log_num + std::log(denom_sum);

      for (int k = 0; k < Ki; ++k) {
        prob(p, idx + k) = std::exp(log_num[k] - log_denom);
      }
    }
    idx += Ki;
  }

  return prob;
}

// [[Rcpp::export]]
NumericMatrix cpp_model_MGPCM(NumericMatrix theta, NumericMatrix par) {
  const int N = theta.nrow();
  const int D = theta.ncol();
  const int I = par.nrow();

  if (D <= 0) {
    stop("'theta' must have at least one column.");
  }
  if (par.ncol() <= D) {
    stop("'par' must have D discrimination columns followed by category intercept columns.");
  }

  std::vector<int> length_poly = get_length_poly_mgpcm(par, D);
  int total_cols = 0;
  for (int Ki : length_poly) {
    total_cols += Ki;
  }

  NumericMatrix prob(N, total_cols);
  int idx = 0;
  std::vector<double> log_num;

  for (int i = 0; i < I; ++i) {
    const int Ki = length_poly[i];
    log_num.assign(Ki, 0.0);

    for (int p = 0; p < N; ++p) {
      double eta = 0.0;
      for (int d = 0; d < D; ++d) {
        eta += theta(p, d) * par(i, d);
      }

      double max_log_num = -std::numeric_limits<double>::infinity();
      for (int k = 0; k < Ki; ++k) {
        log_num[k] = static_cast<double>(k) * eta + par(i, D + k);
        if (log_num[k] > max_log_num) {
          max_log_num = log_num[k];
        }
      }

      double denom_sum = 0.0;
      for (int k = 0; k < Ki; ++k) {
        denom_sum += std::exp(log_num[k] - max_log_num);
      }
      const double log_denom = max_log_num + std::log(denom_sum);

      for (int k = 0; k < Ki; ++k) {
        prob(p, idx + k) = std::exp(log_num[k] - log_denom);
      }
    }

    idx += Ki;
  }

  return prob;
}

// [[Rcpp::export]]
NumericMatrix cpp_model_FCGGUM(NumericMatrix theta, NumericMatrix par,
                               List patterns_total, List patterns) {
  const int D = theta.ncol();
  if (D <= 0) {
    stop("'theta' must have at least one column.");
  }
  NumericMatrix agree = mggum_binary_agree(theta, par);
  return forced_choice_from_agree(agree, patterns_total, patterns);
}

// [[Rcpp::export]]
NumericMatrix cpp_model_FCMIRT(NumericMatrix theta, NumericMatrix par,
                               List patterns_total, List patterns) {
  NumericMatrix agree = mirt_binary_agree(theta, par);
  return forced_choice_from_agree(agree, patterns_total, patterns);
}

// [[Rcpp::export]]
NumericMatrix cpp_model_FCDCM(IntegerMatrix zeta, NumericMatrix par,
                              IntegerMatrix patterns) {
  const int N = zeta.nrow();
  const int I_states = zeta.ncol();
  const int B = patterns.nrow();
  if (patterns.ncol() != 2) {
    stop("'patterns' must have two columns for FCDCM blocks.");
  }
  if (par.nrow() != B || (par.ncol() != 2 && par.ncol() != 4)) {
    stop("'par' must be a B x 2 or B x 4 FCDCM parameter matrix.");
  }

  NumericMatrix prob(N, B);
  std::vector<int> item_a(B), item_b(B);
  for (int b = 0; b < B; ++b) {
    item_a[b] = fcdcm_item_index(patterns(b, 0), I_states);
    item_b[b] = fcdcm_item_index(patterns(b, 1), I_states);
  }

  for (int n = 0; n < N; ++n) {
    for (int b = 0; b < B; ++b) {
      prob(n, b) = fcdcm_block_prob(
        zeta(n, item_a[b]), zeta(n, item_b[b]),
        par(b, 0), par(b, 1),
        par.ncol() == 4 ? par(b, 2) : 0.5,
        par.ncol() == 4 ? par(b, 3) : 0.5
      );
    }
  }
  return prob;
}

// [[Rcpp::export]]
NumericMatrix cpp_model_FCDCM_marginal(NumericMatrix theta,
                                       NumericVector delta1,
                                       NumericVector delta0,
                                       NumericMatrix par,
                                       NumericMatrix alpha_patterns,
                                       IntegerMatrix zeta_patterns,
                                       IntegerMatrix patterns) {
  fcdcm_check_dims(theta, delta1, delta0, par,
                   alpha_patterns, zeta_patterns, patterns);
  const int N = theta.nrow();
  const int C = alpha_patterns.nrow();
  const int B = patterns.nrow();

  std::vector<int> block_case = fcdcm_block_case(zeta_patterns, patterns);
  NumericMatrix class_prob(C, B);
  fcdcm_class_prob(par, block_case, C, B, class_prob);

  NumericMatrix prob(N, B);
  std::vector<double> weight(C);
  for (int n = 0; n < N; ++n) {
    fcdcm_alpha_weights(theta(n, 0), delta1, delta0,
                        alpha_patterns, weight);
    for (int b = 0; b < B; ++b) {
      double pb = 0.0;
      for (int c = 0; c < C; ++c) {
        pb += weight[c] * class_prob(c, b);
      }
      prob(n, b) = fcdcm_clip_prob(pb);
    }
  }
  return prob;
}

// [[Rcpp::export]]
NumericVector cpp_fcdcm_case_loglik_class_prob(NumericMatrix theta,
                                               IntegerMatrix response,
                                               NumericVector delta1,
                                               NumericVector delta0,
                                               NumericMatrix alpha_patterns,
                                               NumericMatrix class_prob) {
  const int N = theta.nrow();
  const int C = alpha_patterns.nrow();
  const int D = alpha_patterns.ncol();
  const int B = class_prob.ncol();
  if (theta.ncol() != 1) {
    stop("FCDCM theta must be an N x 1 matrix.");
  }
  if (delta1.size() != D || delta0.size() != D) {
    stop("'delta1' and 'delta0' must have length equal to ncol(alpha_patterns).");
  }
  if (class_prob.nrow() != C) {
    stop("'class_prob' must have one row per alpha pattern.");
  }
  if (response.nrow() != N || response.ncol() != B) {
    stop("'response' dimensions must match nrow(theta) and ncol(class_prob).");
  }

  std::vector<double> log_p1, log_p0;
  fcdcm_class_logprob_from_prob(class_prob, log_p1, log_p0);

  NumericVector out(N);
  for (int n = 0; n < N; ++n) {
    out[n] = fcdcm_person_loglik(theta(n, 0), response, n,
                                 delta1, delta0, alpha_patterns,
                                 log_p1, log_p0);
  }
  return out;
}

// [[Rcpp::export]]
List cpp_gibbs_fcdcm_theta(NumericMatrix theta,
                            IntegerMatrix response,
                            NumericVector delta1,
                            NumericVector delta0,
                            NumericMatrix par,
                            NumericMatrix alpha_patterns,
                            IntegerMatrix zeta_patterns,
                            IntegerMatrix patterns,
                            double step,
                            double lower,
                            double upper) {
  fcdcm_check_dims(theta, delta1, delta0, par,
                   alpha_patterns, zeta_patterns, patterns);
  const int N = theta.nrow();
  const int C = alpha_patterns.nrow();
  const int B = patterns.nrow();
  if (theta.ncol() != 1) {
    stop("FCDCM theta must be an N x 1 matrix.");
  }
  if (response.nrow() != N || response.ncol() != B) {
    stop("'response' dimensions must match nrow(theta) and nrow(patterns).");
  }

  std::vector<int> block_case = fcdcm_block_case(zeta_patterns, patterns);
  std::vector<double> log_p1, log_p0;
  fcdcm_class_logprob(par, block_case, C, B, log_p1, log_p0);

  const int L = forcechoice_istem::grid_length_from_step(step);
  std::vector<double> axis = forcechoice_istem::make_axis(lower, upper, L);
  std::vector<double> person_class_ll(static_cast<std::size_t>(N) * C);
  for (int n = 0; n < N; ++n) {
    for (int c = 0; c < C; ++c) {
      double val = 0.0;
      for (int b = 0; b < B; ++b) {
        const int y = response(n, b);
        val += (y == 1) ? log_p1[c * B + b] : log_p0[c * B + b];
      }
      person_class_ll[static_cast<std::size_t>(n) * C + c] = val;
    }
  }

  std::vector<double> grid_class_lp(static_cast<std::size_t>(L) * C);
  for (int g = 0; g < L; ++g) {
    for (int c = 0; c < C; ++c) {
      grid_class_lp[static_cast<std::size_t>(g) * C + c] =
        fcdcm_log_alpha_prob(axis[g], alpha_patterns, c, delta1, delta0);
    }
  }

  std::vector<double> logw(L);
  std::vector<double> class_lp(C);

  NumericMatrix out = clone(theta);
  int sampled = 0;
  for (int n = 0; n < N; ++n) {
    for (int g = 0; g < L; ++g) {
      double value = axis[g];
      for (int c = 0; c < C; ++c) {
        class_lp[c] = grid_class_lp[static_cast<std::size_t>(g) * C + c] +
          person_class_ll[static_cast<std::size_t>(n) * C + c];
      }
      logw[g] = logspace_sum(class_lp) - 0.5 * value * value;
    }
    int idx = forcechoice_istem::sample_log_weight(logw);
    if (idx >= 0) {
      out(n, 0) = axis[idx];
      ++sampled;
    }
  }

  return List::create(
    _["theta"] = out,
    _["accepted"] = sampled,
    _["proposed"] = N
  );
}

// [[Rcpp::export]]
List cpp_fcdcm_sample_alpha(NumericMatrix theta,
                            IntegerMatrix response,
                            NumericVector delta1,
                            NumericVector delta0,
                            NumericMatrix alpha_patterns,
                            IntegerMatrix zeta_patterns,
                            NumericMatrix class_prob) {
  const int N = theta.nrow();
  const int C = alpha_patterns.nrow();
  const int D = alpha_patterns.ncol();
  const int I = zeta_patterns.ncol();
  const int B = class_prob.ncol();
  if (theta.ncol() != 1) {
    stop("FCDCM theta must be an N x 1 matrix.");
  }
  if (delta1.size() != D || delta0.size() != D) {
    stop("'delta1' and 'delta0' must have length equal to ncol(alpha_patterns).");
  }
  if (zeta_patterns.nrow() != C || class_prob.nrow() != C) {
    stop("'zeta_patterns' and 'class_prob' must have one row per alpha pattern.");
  }
  if (response.nrow() != N || response.ncol() != B) {
    stop("'response' dimensions must match nrow(theta) and ncol(class_prob).");
  }

  std::vector<double> log_p1, log_p0;
  fcdcm_class_logprob_from_prob(class_prob, log_p1, log_p0);

  IntegerMatrix alpha(N, D);
  IntegerMatrix zeta(N, I);
  IntegerVector class_index(N);
  std::vector<double> lp(C);

  for (int n = 0; n < N; ++n) {
    for (int c = 0; c < C; ++c) {
      double val = fcdcm_log_alpha_prob(theta(n, 0), alpha_patterns, c,
                                        delta1, delta0);
      for (int b = 0; b < B; ++b) {
        const int y = response(n, b);
        val += (y == 1) ? log_p1[c * B + b] : log_p0[c * B + b];
      }
      lp[c] = val;
    }

    const double log_denom = logspace_sum(lp);
    double u = R::runif(0.0, 1.0);
    double cum = 0.0;
    int chosen = C - 1;
    for (int c = 0; c < C; ++c) {
      cum += std::exp(lp[c] - log_denom);
      if (u <= cum) {
        chosen = c;
        break;
      }
    }

    class_index[n] = chosen + 1;
    for (int d = 0; d < D; ++d) {
      alpha(n, d) = (alpha_patterns(chosen, d) > 0.5) ? 1 : 0;
    }
    for (int i = 0; i < I; ++i) {
      zeta(n, i) = zeta_patterns(chosen, i);
    }
  }

  return List::create(
    _["alpha"] = alpha,
    _["zeta"] = zeta,
    _["class.index"] = class_index
  );
}

NumericMatrix cpp_fcdcm_posterior_alpha(NumericVector theta,
                                        IntegerMatrix response,
                                        NumericVector delta1,
                                        NumericVector delta0,
                                        NumericMatrix par,
                                        NumericMatrix alpha_patterns,
                                        IntegerMatrix zeta_patterns,
                                        IntegerMatrix patterns) {
  NumericMatrix theta_mat(theta.size(), 1);
  for (int n = 0; n < theta.size(); ++n) theta_mat(n, 0) = theta[n];
  fcdcm_check_dims(theta_mat, delta1, delta0, par,
                   alpha_patterns, zeta_patterns, patterns);
  const int N = theta.size();
  const int C = alpha_patterns.nrow();
  const int D = alpha_patterns.ncol();
  const int B = patterns.nrow();
  if (response.nrow() != N || response.ncol() != B) {
    stop("'response' dimensions must match theta and nrow(patterns).");
  }

  std::vector<int> block_case = fcdcm_block_case(zeta_patterns, patterns);
  std::vector<double> log_p1, log_p0;
  fcdcm_class_logprob(par, block_case, C, B, log_p1, log_p0);

  NumericMatrix out(N, D);
  std::vector<double> lp(C);
  for (int n = 0; n < N; ++n) {
    for (int c = 0; c < C; ++c) {
      double val = fcdcm_log_alpha_prob(theta[n], alpha_patterns, c,
                                        delta1, delta0);
      for (int b = 0; b < B; ++b) {
        const int y = response(n, b);
        val += (y == 1) ? log_p1[c * B + b] : log_p0[c * B + b];
      }
      lp[c] = val;
    }
    const double log_denom = logspace_sum(lp);
    for (int c = 0; c < C; ++c) {
      const double w = std::exp(lp[c] - log_denom);
      for (int d = 0; d < D; ++d) {
        out(n, d) += w * alpha_patterns(c, d);
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix cpp_fcdcm_posterior_alpha_class_prob(NumericVector theta,
                                                   IntegerMatrix response,
                                                   NumericVector delta1,
                                                   NumericVector delta0,
                                                   NumericMatrix alpha_patterns,
                                                   NumericMatrix class_prob) {
  const int N = theta.size();
  const int C = alpha_patterns.nrow();
  const int D = alpha_patterns.ncol();
  const int B = class_prob.ncol();
  if (delta1.size() != D || delta0.size() != D) {
    stop("'delta1' and 'delta0' must have length equal to ncol(alpha_patterns).");
  }
  if (class_prob.nrow() != C) {
    stop("'class_prob' must have one row per alpha pattern.");
  }
  if (response.nrow() != N || response.ncol() != B) {
    stop("'response' dimensions must match theta and ncol(class_prob).");
  }

  std::vector<double> log_p1, log_p0;
  fcdcm_class_logprob_from_prob(class_prob, log_p1, log_p0);

  NumericMatrix out(N, D);
  std::vector<double> lp(C);
  for (int n = 0; n < N; ++n) {
    for (int c = 0; c < C; ++c) {
      double val = fcdcm_log_alpha_prob(theta[n], alpha_patterns, c,
                                        delta1, delta0);
      for (int b = 0; b < B; ++b) {
        const int y = response(n, b);
        val += (y == 1) ? log_p1[c * B + b] : log_p0[c * B + b];
      }
      lp[c] = val;
    }
    const double log_denom = logspace_sum(lp);
    for (int c = 0; c < C; ++c) {
      const double w = std::exp(lp[c] - log_denom);
      for (int d = 0; d < D; ++d) {
        out(n, d) += w * alpha_patterns(c, d);
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
List cpp_fcdcm_latent_support(NumericMatrix theta,
                              NumericVector pi,
                              NumericVector delta1,
                              NumericVector delta0,
                              NumericMatrix par,
                              NumericMatrix alpha_patterns,
                              IntegerMatrix zeta_patterns,
                              IntegerMatrix patterns) {
  fcdcm_check_dims(theta, delta1, delta0, par,
                   alpha_patterns, zeta_patterns, patterns);
  const int L = theta.nrow();
  const int C = alpha_patterns.nrow();
  const int D = alpha_patterns.ncol();
  const int B = patterns.nrow();
  if (pi.size() != L) {
    stop("'pi' must have one weight per theta support point.");
  }

  double pi_sum = 0.0;
  for (int l = 0; l < L; ++l) {
    if (pi[l] > 0.0) pi_sum += pi[l];
  }
  if (pi_sum <= 0.0 || !std::isfinite(pi_sum)) {
    stop("'pi' must contain positive finite weights.");
  }

  std::vector<int> block_case = fcdcm_block_case(zeta_patterns, patterns);
  NumericMatrix class_prob(C, B);
  fcdcm_class_prob(par, block_case, C, B, class_prob);

  NumericMatrix prob(L * C, B);
  NumericVector support_pi(L * C);
  NumericMatrix support_theta(L * C, 1);
  NumericMatrix support_alpha(L * C, D);
  NumericMatrix alpha_weight(L, C);
  std::vector<double> weight(C);

  for (int l = 0; l < L; ++l) {
    fcdcm_alpha_weights(theta(l, 0), delta1, delta0,
                        alpha_patterns, weight);
    const double pi_l = std::max(pi[l], 0.0) / pi_sum;
    for (int c = 0; c < C; ++c) {
      const int row = l * C + c;
      alpha_weight(l, c) = weight[c];
      support_pi[row] = pi_l * weight[c];
      support_theta(row, 0) = theta(l, 0);
      for (int d = 0; d < D; ++d) {
        support_alpha(row, d) = alpha_patterns(c, d);
      }
      for (int b = 0; b < B; ++b) {
        prob(row, b) = class_prob(c, b);
      }
    }
  }

  return List::create(
    _["prob"] = prob,
    _["pi"] = support_pi,
    _["theta"] = support_theta,
    _["alpha"] = support_alpha,
    _["alpha.weight"] = alpha_weight,
    _["class.prob"] = class_prob
  );
}

// ---- FCGDINA model -----------------------------------------------------------

namespace {

// Compute per-class item endorsement probability matrix.
// alpha_pat: C x D matrix of attribute patterns
// flat_design: flattened design matrices (row-major per item)
// design_rows, design_cols, design_offset: per-item design matrix specs
// delta: stacked delta vector (same offset scheme)
// class_map: for each item, which row of design matrix each class maps to
NumericMatrix fcgdina_item_prob(const IntegerMatrix& alpha_pat,
                                const NumericVector& flat_design,
                                const IntegerVector& design_rows,
                                const IntegerVector& design_cols,
                                const IntegerVector& design_offset,
                                const IntegerVector& delta_offset,
                                const NumericVector& delta,
                                const List& class_map_list,
                                int C, int I_states,
                                double eps = 1e-12) {
  NumericMatrix prob(C, I_states);

  for (int i = 0; i < I_states; ++i) {
    int off_x = design_offset[i];  // offset in flat_design entries
    int off_d = delta_offset[i];   // offset in delta entries
    int nr = design_rows[i];
    int nc = design_cols[i];
    if (nr <= 0 || nc <= 0) continue;

    IntegerVector map_i = as<IntegerVector>(class_map_list[i]);

    // Pre-compute raw item probabilities: X * delta
    std::vector<double> p_raw(nr);
    for (int r = 0; r < nr; ++r) {
      double val = 0.0;
      for (int c = 0; c < nc; ++c) {
        val += flat_design[off_x + r * nc + c] * delta[off_d + c];
      }
      // Clamp
      if (val < eps) val = eps;
      if (val > 1.0 - eps) val = 1.0 - eps;
      p_raw[r] = val;
    }

    // Map to classes
    for (int cls = 0; cls < C; ++cls) {
      int row = map_i[cls] - 1;  // 1-based to 0-based
      if (row < 0) row = 0;
      if (row >= nr) row = nr - 1;
      prob(cls, i) = p_raw[row];
    }
  }

  return prob;
}

} // namespace

// [[Rcpp::export]]
NumericMatrix cpp_model_FCGDINA(const IntegerMatrix& alpha_pat,
                                const NumericVector& flat_design,
                                const IntegerVector& design_rows,
                                const IntegerVector& design_cols,
                                const IntegerVector& design_offset,
                                const IntegerVector& delta_offset,
                                const NumericVector& delta,
                                const List& class_map_list,
                                const List& patterns_total,
                                const List& patterns,
                                int C, int I_states) {
  NumericMatrix agree = fcgdina_item_prob(
    alpha_pat, flat_design, design_rows, design_cols,
    design_offset, delta_offset,
    delta, class_map_list, C, I_states
  );
  return forced_choice_from_agree(agree, patterns_total, patterns);
}

// [[Rcpp::export]]
NumericMatrix cpp_fcgdina_class_posterior(const NumericMatrix& prob_class,
                                          const IntegerMatrix& response,
                                          const List& patterns,
                                          const NumericVector& prior,
                                          int N, int B, int C) {
  // prob_class: C x sum(pattern_counts) matrix
  // response: N x B integer matrix
  // prior: length-C vector

  NumericMatrix post(N, C);

  // Compute pattern column offsets per block
  int n_blocks = patterns.size();
  std::vector<int> block_offset(n_blocks + 1, 0);
  for (int b = 0; b < n_blocks; ++b) {
    NumericMatrix pat = as<NumericMatrix>(patterns[b]);
    block_offset[b + 1] = block_offset[b] + pat.nrow();
  }

  for (int n = 0; n < N; ++n) {
    std::vector<double> lp(C, 0.0);
    int col_idx = 0;

    for (int b = 0; b < B; ++b) {
      NumericMatrix pat = as<NumericMatrix>(patterns[b]);
      int n_pat = pat.nrow();
      int y = response(n, b) - 1;  // 0-based

      if (y >= 0 && y < n_pat) {
        for (int c = 0; c < C; ++c) {
          double p = prob_class(c, col_idx + y);
          if (p < 1e-16) p = 1e-16;
          lp[c] += std::log(p);
        }
      }
      col_idx += n_pat;
    }

    // Add log prior and normalize
    double max_lp = -std::numeric_limits<double>::infinity();
    for (int c = 0; c < C; ++c) {
      lp[c] += std::log(prior[c] > 1e-16 ? prior[c] : 1e-16);
      if (lp[c] > max_lp) max_lp = lp[c];
    }

    double denom = 0.0;
    for (int c = 0; c < C; ++c) {
      double v = std::exp(lp[c] - max_lp);
      post(n, c) = v;
      denom += v;
    }
    for (int c = 0; c < C; ++c) {
      post(n, c) /= denom;
    }
  }

  return post;
}

// [[Rcpp::export]]
double cpp_fcgdina_weighted_loglik(const NumericMatrix& prob_class,
                                   const IntegerMatrix& response,
                                   const List& patterns,
                                   const NumericMatrix& class_weight,
                                   int N, int B, int C) {
  // prob_class: C x sum(pattern_counts)
  // class_weight: N x C posterior weights or stochastic class indicators
  double ell = 0.0;
  int col_idx = 0;

  for (int b = 0; b < B; ++b) {
    NumericMatrix pat = as<NumericMatrix>(patterns[b]);
    int n_pat = pat.nrow();

    for (int n = 0; n < N; ++n) {
      int y = response(n, b) - 1;
      if (y < 0 || y >= n_pat) continue;

      for (int c = 0; c < C; ++c) {
        double w = class_weight(n, c);
        if (w == 0.0) continue;
        double p = prob_class(c, col_idx + y);
        if (p < 1e-16) p = 1e-16;
        ell += w * std::log(p);
      }
    }

    col_idx += n_pat;
  }

  return ell;
}

// [[Rcpp::export]]
double cpp_fcgdina_weighted_loglik_block(const NumericMatrix& prob_block,
                                         const IntegerVector& response_block,
                                         const NumericMatrix& class_weight,
                                         int N, int C) {
  // prob_block: C x observed-pattern-count for one block
  double ell = 0.0;
  const int n_pat = prob_block.ncol();

  for (int n = 0; n < N; ++n) {
    int y = response_block[n] - 1;
    if (y < 0 || y >= n_pat) continue;

    for (int c = 0; c < C; ++c) {
      double w = class_weight(n, c);
      if (w == 0.0) continue;
      double p = prob_block(c, y);
      if (p < 1e-16) p = 1e-16;
      ell += w * std::log(p);
    }
  }

  return ell;
}

// [[Rcpp::export]]
double cpp_fcgdina_marginal_loglik(const NumericMatrix& prob_class,
                                   const IntegerMatrix& response,
                                   const List& patterns,
                                   const NumericVector& prior,
                                   int N, int B, int C) {
  double lik = 0.0;

  for (int n = 0; n < N; ++n) {
    std::vector<double> lp(C, 0.0);
    int col_idx = 0;

    for (int b = 0; b < B; ++b) {
      NumericMatrix pat = as<NumericMatrix>(patterns[b]);
      int n_pat = pat.nrow();
      int y = response(n, b) - 1;

      if (y >= 0 && y < n_pat) {
        for (int c = 0; c < C; ++c) {
          double p = prob_class(c, col_idx + y);
          if (p < 1e-16) p = 1e-16;
          lp[c] += std::log(p);
        }
      }

      col_idx += n_pat;
    }

    double max_lp = -std::numeric_limits<double>::infinity();
    for (int c = 0; c < C; ++c) {
      lp[c] += std::log(prior[c] > 1e-16 ? prior[c] : 1e-16);
      if (lp[c] > max_lp) max_lp = lp[c];
    }

    double denom = 0.0;
    for (int c = 0; c < C; ++c) {
      denom += std::exp(lp[c] - max_lp);
    }
    lik += max_lp + std::log(denom);
  }

  return lik;
}

// ---- FCGDINA logit-to-probability delta conversion -------------------------
// Converts delta from logit scale (Stan parameterisation) to probability
// scale (EM / iStEM parameterisation).  The transformation is applied to
// every posterior draw so that posterior summaries (mean, SE, Rhat) are
// reported on the natural probability scale.
//
// For item i with design matrix X_i (nr × nc), logit-scale delta δ_logit,
// and endorsement probabilities p = σ(X_i · δ_logit), the probability-scale
// delta is the least-squares solution to  X_i · δ_prob ≈ p:
//
//   δ_prob  =  (X_i' X_i)⁻¹  X_i'  σ(X_i · δ_logit)
//
// When X_i is square (GDINA) this is the exact inverse; for reduced-rank
// designs (DINA, DINO, ACDM) it gives the constrained CDM parameters that
// best reproduce the endorsement probabilities.
//
// The function loops over draws externally so that the per-item design
// matrices are unpacked only once.

namespace {

// Convert one draw.  delta_prob is already sized and we overwrite in-place.
void fcgdina_delta_logit_to_prob_one(
    const double* delta_logit,   // length total_delta_len
    double* delta_prob,          // length total_delta_len (output)
    const double* flat_design,
    const int* design_rows,
    const int* design_cols,
    const int* design_offset,
    const int* delta_offset,
    int I_states) {

  for (int i = 0; i < I_states; ++i) {
    int off_x = design_offset[i];
    int off_d = delta_offset[i];
    int nr = design_rows[i];
    int nc = design_cols[i];

    if (nr <= 0 || nc <= 0) continue;

    // Build X (nr × nc, row-major) and compute eta = X * delta_logit
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                   Eigen::RowMajor>>
        X(flat_design + off_x, nr, nc);

    Eigen::Map<const Eigen::VectorXd> dlogit(delta_logit + off_d, nc);
    Eigen::VectorXd eta = X * dlogit;

    // p = σ(eta), clamped
    Eigen::VectorXd p(nr);
    for (int r = 0; r < nr; ++r) {
      double x = eta(r);
      double pr;
      if (x >= 0.0) {
        pr = 1.0 / (1.0 + std::exp(-x));
      } else {
        double z = std::exp(x);
        pr = z / (1.0 + z);
      }
      if (pr < 1e-12) pr = 1e-12;
      if (pr > 1.0 - 1e-12) pr = 1.0 - 1e-12;
      p(r) = pr;
    }

    // Least squares: solve  (X'X) · dprob = X'p
    Eigen::MatrixXd XtX = X.transpose() * X;
    Eigen::VectorXd Xtp = X.transpose() * p;
    Eigen::VectorXd dprob = XtX.ldlt().solve(Xtp);

    for (int c = 0; c < nc; ++c) {
      delta_prob[off_d + c] = dprob(c);
    }
  }
}

} // namespace

// [[Rcpp::export]]
NumericMatrix cpp_fcgdina_delta_logit_to_prob(
    const NumericMatrix& delta_logit,
    const NumericVector& flat_design,
    const IntegerVector& design_rows,
    const IntegerVector& design_cols,
    const IntegerVector& design_offset,
    const IntegerVector& delta_offset,
    int I_states) {

  int n_draws = delta_logit.nrow();
  int total_len = delta_logit.ncol();
  NumericMatrix out(n_draws, total_len);

  // Extract raw pointers for speed
  const double* fd = REAL(flat_design);
  const int* dr = INTEGER(design_rows);
  const int* dc = INTEGER(design_cols);
  const int* dx = INTEGER(design_offset);
  const int* dd = INTEGER(delta_offset);

  std::vector<double> src(total_len);
  std::vector<double> dst(total_len);
  for (int draw = 0; draw < n_draws; ++draw) {
    for (int j = 0; j < total_len; ++j) {
      src[j] = delta_logit(draw, j);
    }
    fcgdina_delta_logit_to_prob_one(
        src.data(), dst.data(), fd, dr, dc, dx, dd, I_states);
    for (int j = 0; j < total_len; ++j) {
      out(draw, j) = dst[j];
    }
  }

  return out;
}
