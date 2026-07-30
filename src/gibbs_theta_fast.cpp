// Block Gibbs theta samplers on a finite D-dimensional support.
// Each respondent is sampled from p(theta_i | all responses_i, item, Corr)
// by normalizing that respondent's complete log posterior over the grid.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "istem_gibbs_grid.h"
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

constexpr double kProbFloor = 1e-12;
constexpr double kMirtFloor = 1e-10;

using forcechoice_istem::ThetaGrid;
using forcechoice_istem::build_theta_grid;
using forcechoice_istem::logspace_add;
using forcechoice_istem::logspace_sum;
using forcechoice_istem::sample_log_weight;

inline double inv_logit(double x) {
  if (x >= 0.0) {
    double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  double z = std::exp(x);
  return z / (1.0 + z);
}

template <typename MatrixType>
void set_theta_row(MatrixType& theta, int person, const ThetaGrid& grid,
                   int grid_index) {
  const double* x = &grid.point[static_cast<size_t>(grid_index) * grid.D];
  for (int d = 0; d < grid.D; ++d) theta(person, d) = x[d];
}

double mirt_person_loglik(const double* x, int person, const MatrixXd& par,
                          const IntegerMatrix& response, int D, int I) {
  double out = 0.0;
  for (int j = 0; j < I; ++j) {
    double eta = -par(j, D);
    for (int d = 0; d < D; ++d) eta += x[d] * par(j, d);
    double p = inv_logit(eta);
    double c = par(j, D + 1);
    double dval = par(j, D + 2);
    p = c + (dval - c) * p;
    if (p < kMirtFloor) p = kMirtFloor;
    if (p > 1.0 - kMirtFloor) p = 1.0 - kMirtFloor;
    int y = response(person, j);
    out += (y == 1) ? std::log(p) : std::log1p(-p);
  }
  return out;
}

double mgpcm_person_loglik(const double* x, int person, const MatrixXd& par,
                           const IntegerMatrix& response,
                           const IntegerVector& length_poly, int D, int I) {
  double out = 0.0;
  for (int j = 0; j < I; ++j) {
    double eta = 0.0;
    for (int d = 0; d < D; ++d) eta += x[d] * par(j, d);
    int Ki = length_poly[j];
    int y = response(person, j);
    if (y < 0 || y >= Ki) continue;

    std::vector<double> lp(Ki);
    for (int k = 0; k < Ki; ++k) lp[k] = k * eta + par(j, D + k);
    out += lp[y] - logspace_sum(lp);
  }
  return out;
}

double mggum_person_loglik(const double* x, int person, const MatrixXd& par,
                           const IntegerMatrix& response,
                           const IntegerVector& length_poly, int D, int I) {
  double out = 0.0;
  for (int j = 0; j < I; ++j) {
    int Ki = length_poly[j];
    int y = response(person, j);
    if (y < 0 || y >= Ki) continue;

    double sum_a = 0.0;
    double rsq = 0.0;
    for (int d = 0; d < D; ++d) {
      sum_a += par(j, d);
      double diff = x[d] - par(j, D + d);
      rsq += par(j, d) * par(j, d) * diff * diff;
    }
    double r = std::sqrt(rsq + 1e-10);
    int Mi = 2 * Ki - 1;
    double cp = 0.0;
    std::vector<double> lp(Ki);
    for (int k = 0; k < Ki; ++k) {
      cp += par(j, 2 * D + k) * sum_a;
      double zk = static_cast<double>(k);
      lp[k] = logspace_add(zk * r - cp,
                           (static_cast<double>(Mi) - zk) * r - cp);
    }
    out += lp[y] - logspace_sum(lp);
  }
  return out;
}

struct FcBlock {
  int tr;
  int tc;
  int orr;
  int oc;
  std::vector<int> total;
  std::vector<int> observed;
  std::vector<std::vector<int>> map;
};

FcBlock parse_fc(SEXP patterns_total, SEXP patterns) {
  NumericMatrix total(patterns_total);
  NumericMatrix observed(patterns);
  FcBlock b;
  b.tr = total.nrow();
  b.tc = total.ncol();
  b.orr = observed.nrow();
  b.oc = observed.ncol();
  b.total.resize(b.tr * b.tc);
  b.observed.resize(b.orr * b.oc);

  for (int r = 0; r < b.tr; ++r) {
    for (int c = 0; c < b.tc; ++c) b.total[r * b.tc + c] = total(r, c) - 1;
  }
  for (int r = 0; r < b.orr; ++r) {
    for (int c = 0; c < b.oc; ++c) b.observed[r * b.oc + c] = observed(r, c) - 1;
  }

  b.map.resize(b.orr);
  if (b.oc == b.tc) {
    for (int o = 0; o < b.orr; ++o) {
      for (int t = 0; t < b.tr; ++t) {
        bool ok = true;
        for (int c = 0; c < b.tc; ++c) {
          if (b.observed[o * b.oc + c] != b.total[t * b.tc + c]) {
            ok = false;
            break;
          }
        }
        if (ok) b.map[o].push_back(t);
      }
    }
  } else if (b.oc == 2) {
    for (int o = 0; o < b.orr; ++o) {
      for (int t = 0; t < b.tr; ++t) {
        if (b.observed[o * b.oc] == b.total[t * b.tc] &&
            b.observed[o * b.oc + 1] == b.total[t * b.tc + b.tc - 1]) {
          b.map[o].push_back(t);
        }
      }
    }
  } else if (b.oc == 1) {
    for (int o = 0; o < b.orr; ++o) {
      for (int t = 0; t < b.tr; ++t) {
        if (b.observed[o * b.oc] == b.total[t * b.tc]) b.map[o].push_back(t);
      }
    }
  }
  return b;
}

double fc_block_loglik(const std::vector<double>& item_logit,
                       const FcBlock& b, int response_index) {
  std::vector<double> full_lp(b.tr);
  for (int r = 0; r < b.tr; ++r) {
    double lp = 0.0;
    for (int pos = 0; pos < b.tc - 1; ++pos) {
      int selected = b.total[r * b.tc + pos];
      double mx = -std::numeric_limits<double>::infinity();
      for (int cp = pos; cp < b.tc; ++cp) {
        int item = b.total[r * b.tc + cp];
        if (item_logit[item] > mx) mx = item_logit[item];
      }
      double denom = 0.0;
      for (int cp = pos; cp < b.tc; ++cp) {
        int item = b.total[r * b.tc + cp];
        denom += std::exp(item_logit[item] - mx);
      }
      lp += item_logit[selected] - (mx + std::log(denom));
    }
    full_lp[r] = lp;
  }

  std::vector<double> observed_lp(b.orr);
  for (int o = 0; o < b.orr; ++o) {
    double mx = -std::numeric_limits<double>::infinity();
    for (int total_index : b.map[o]) {
      if (full_lp[total_index] > mx) mx = full_lp[total_index];
    }
    double s = 0.0;
    for (int total_index : b.map[o]) s += std::exp(full_lp[total_index] - mx);
    observed_lp[o] = mx + std::log(s + 1e-15);
  }
  return observed_lp[response_index] - logspace_sum(observed_lp);
}

double fcmirt_person_loglik(const double* x, int person, const MatrixXd& par,
                            const IntegerMatrix& response,
                            const std::vector<FcBlock>& blocks,
                            int D, int I, int n_blocks) {
  std::vector<double> item_logit(I);
  for (int j = 0; j < I; ++j) {
    double eta = -par(j, D);
    for (int d = 0; d < D; ++d) eta += x[d] * par(j, d);
    double p = inv_logit(eta);
    double c = par(j, D + 1);
    double dval = par(j, D + 2);
    p = c + (dval - c) * p;
    if (p < kProbFloor) p = kProbFloor;
    if (p > 1.0 - kProbFloor) p = 1.0 - kProbFloor;
    item_logit[j] = std::log(p) - std::log1p(-p);
  }

  double out = 0.0;
  for (int b = 0; b < n_blocks; ++b) {
    int y = response(person, b);
    if (y >= 1) out += fc_block_loglik(item_logit, blocks[b], y - 1);
  }
  return out;
}

double fcggum_person_loglik(const double* x, int person, const MatrixXd& par,
                            const IntegerMatrix& response,
                            const std::vector<FcBlock>& blocks,
                            int D, int I, int n_blocks) {
  std::vector<double> item_logit(I);
  for (int j = 0; j < I; ++j) {
    double sum_a = 0.0;
    double rsq = 0.0;
    for (int d = 0; d < D; ++d) {
      sum_a += par(j, d);
      double diff = x[d] - par(j, D + d);
      rsq += par(j, d) * par(j, d) * diff * diff;
    }
    double r = std::sqrt(rsq + 1e-10);
    double c0 = par(j, 2 * D) * sum_a;
    double c1 = c0 + par(j, 2 * D + 1) * sum_a;
    double n0 = logspace_add(-c0, 3 * r - c0);
    double n1 = logspace_add(r - c1, 2 * r - c1);
    double lv = n1 - n0;
    if (lv > 27.63) lv = 27.63;
    if (lv < -27.63) lv = -27.63;
    item_logit[j] = lv;
  }

  double out = 0.0;
  for (int b = 0; b < n_blocks; ++b) {
    int y = response(person, b);
    if (y >= 1) out += fc_block_loglik(item_logit, blocks[b], y - 1);
  }
  return out;
}

struct BinaryLogProbTable {
  std::vector<double> logp0;
  std::vector<double> logp1;
};

BinaryLogProbTable mirt_logprob_table(const ThetaGrid& grid,
                                      const MatrixXd& par,
                                      int D, int I) {
  BinaryLogProbTable table;
  table.logp0.resize(static_cast<size_t>(grid.G) * I);
  table.logp1.resize(static_cast<size_t>(grid.G) * I);
  for (int g = 0; g < grid.G; ++g) {
    const double* x = &grid.point[static_cast<size_t>(g) * D];
    const size_t base = static_cast<size_t>(g) * I;
    for (int j = 0; j < I; ++j) {
      double eta = -par(j, D);
      for (int d = 0; d < D; ++d) eta += x[d] * par(j, d);
      double p = inv_logit(eta);
      double c = par(j, D + 1);
      double dval = par(j, D + 2);
      p = c + (dval - c) * p;
      if (p < kMirtFloor) p = kMirtFloor;
      if (p > 1.0 - kMirtFloor) p = 1.0 - kMirtFloor;
      table.logp1[base + j] = std::log(p);
      table.logp0[base + j] = std::log1p(-p);
    }
  }
  return table;
}

std::vector<int> category_offsets(const IntegerVector& length_poly, int I) {
  std::vector<int> offsets(I + 1, 0);
  for (int j = 0; j < I; ++j) {
    offsets[j + 1] = offsets[j] + length_poly[j];
  }
  return offsets;
}

std::vector<double> mgpcm_logprob_table(const ThetaGrid& grid,
                                        const MatrixXd& par,
                                        const IntegerVector& length_poly,
                                        const std::vector<int>& offsets,
                                        int total_cat,
                                        int D, int I) {
  std::vector<double> table(static_cast<size_t>(grid.G) * total_cat);
  for (int g = 0; g < grid.G; ++g) {
    const double* x = &grid.point[static_cast<size_t>(g) * D];
    const size_t base = static_cast<size_t>(g) * total_cat;
    for (int j = 0; j < I; ++j) {
      double eta = 0.0;
      for (int d = 0; d < D; ++d) eta += x[d] * par(j, d);
      int Ki = length_poly[j];
      std::vector<double> lp(Ki);
      for (int k = 0; k < Ki; ++k) lp[k] = k * eta + par(j, D + k);
      double denom = logspace_sum(lp);
      for (int k = 0; k < Ki; ++k) {
        table[base + offsets[j] + k] = lp[k] - denom;
      }
    }
  }
  return table;
}

std::vector<double> mggum_logprob_table(const ThetaGrid& grid,
                                        const MatrixXd& par,
                                        const IntegerVector& length_poly,
                                        const std::vector<int>& offsets,
                                        int total_cat,
                                        int D, int I) {
  std::vector<double> table(static_cast<size_t>(grid.G) * total_cat);
  for (int g = 0; g < grid.G; ++g) {
    const double* x = &grid.point[static_cast<size_t>(g) * D];
    const size_t base = static_cast<size_t>(g) * total_cat;
    for (int j = 0; j < I; ++j) {
      int Ki = length_poly[j];
      double sum_a = 0.0;
      double rsq = 0.0;
      for (int d = 0; d < D; ++d) {
        sum_a += par(j, d);
        double diff = x[d] - par(j, D + d);
        rsq += par(j, d) * par(j, d) * diff * diff;
      }
      double r = std::sqrt(rsq + 1e-10);
      int Mi = 2 * Ki - 1;
      double cp = 0.0;
      std::vector<double> lp(Ki);
      for (int k = 0; k < Ki; ++k) {
        cp += par(j, 2 * D + k) * sum_a;
        double zk = static_cast<double>(k);
        lp[k] = logspace_add(zk * r - cp,
                             (static_cast<double>(Mi) - zk) * r - cp);
      }
      double denom = logspace_sum(lp);
      for (int k = 0; k < Ki; ++k) {
        table[base + offsets[j] + k] = lp[k] - denom;
      }
    }
  }
  return table;
}

std::vector<int> fc_response_offsets(const std::vector<FcBlock>& blocks) {
  const int n_blocks = static_cast<int>(blocks.size());
  std::vector<int> offsets(n_blocks + 1, 0);
  for (int b = 0; b < n_blocks; ++b) {
    offsets[b + 1] = offsets[b] + blocks[b].orr;
  }
  return offsets;
}

void fc_block_logprob(const std::vector<double>& item_logit,
                      const FcBlock& b,
                      std::vector<double>& out) {
  std::vector<double> full_lp(b.tr);
  for (int r = 0; r < b.tr; ++r) {
    double lp = 0.0;
    for (int pos = 0; pos < b.tc - 1; ++pos) {
      int selected = b.total[r * b.tc + pos];
      double mx = -std::numeric_limits<double>::infinity();
      for (int cp = pos; cp < b.tc; ++cp) {
        int item = b.total[r * b.tc + cp];
        if (item_logit[item] > mx) mx = item_logit[item];
      }
      double denom = 0.0;
      for (int cp = pos; cp < b.tc; ++cp) {
        int item = b.total[r * b.tc + cp];
        denom += std::exp(item_logit[item] - mx);
      }
      lp += item_logit[selected] - (mx + std::log(denom));
    }
    full_lp[r] = lp;
  }

  out.assign(b.orr, -std::numeric_limits<double>::infinity());
  for (int o = 0; o < b.orr; ++o) {
    if (b.map[o].empty()) continue;
    double mx = -std::numeric_limits<double>::infinity();
    for (int total_index : b.map[o]) {
      if (full_lp[total_index] > mx) mx = full_lp[total_index];
    }
    double s = 0.0;
    for (int total_index : b.map[o]) s += std::exp(full_lp[total_index] - mx);
    out[o] = mx + std::log(s + 1e-15);
  }
  double denom = logspace_sum(out);
  for (int o = 0; o < b.orr; ++o) out[o] -= denom;
}

std::vector<double> fcmirt_logprob_table(const ThetaGrid& grid,
                                         const MatrixXd& par,
                                         const std::vector<FcBlock>& blocks,
                                         const std::vector<int>& offsets,
                                         int total_resp,
                                         int D, int I) {
  std::vector<double> table(static_cast<size_t>(grid.G) * total_resp);
  std::vector<double> item_logit(I);
  std::vector<double> block_lp;
  for (int g = 0; g < grid.G; ++g) {
    const double* x = &grid.point[static_cast<size_t>(g) * D];
    const size_t base = static_cast<size_t>(g) * total_resp;
    for (int j = 0; j < I; ++j) {
      double eta = -par(j, D);
      for (int d = 0; d < D; ++d) eta += x[d] * par(j, d);
      double p = inv_logit(eta);
      double c = par(j, D + 1);
      double dval = par(j, D + 2);
      p = c + (dval - c) * p;
      if (p < kProbFloor) p = kProbFloor;
      if (p > 1.0 - kProbFloor) p = 1.0 - kProbFloor;
      item_logit[j] = std::log(p) - std::log1p(-p);
    }
    for (int b = 0; b < static_cast<int>(blocks.size()); ++b) {
      fc_block_logprob(item_logit, blocks[b], block_lp);
      for (int o = 0; o < blocks[b].orr; ++o) {
        table[base + offsets[b] + o] = block_lp[o];
      }
    }
  }
  return table;
}

std::vector<double> fcggum_logprob_table(const ThetaGrid& grid,
                                         const MatrixXd& par,
                                         const std::vector<FcBlock>& blocks,
                                         const std::vector<int>& offsets,
                                         int total_resp,
                                         int D, int I) {
  std::vector<double> table(static_cast<size_t>(grid.G) * total_resp);
  std::vector<double> item_logit(I);
  std::vector<double> block_lp;
  for (int g = 0; g < grid.G; ++g) {
    const double* x = &grid.point[static_cast<size_t>(g) * D];
    const size_t base = static_cast<size_t>(g) * total_resp;
    for (int j = 0; j < I; ++j) {
      double sum_a = 0.0;
      double rsq = 0.0;
      for (int d = 0; d < D; ++d) {
        sum_a += par(j, d);
        double diff = x[d] - par(j, D + d);
        rsq += par(j, d) * par(j, d) * diff * diff;
      }
      double r = std::sqrt(rsq + 1e-10);
      double c0 = par(j, 2 * D) * sum_a;
      double c1 = c0 + par(j, 2 * D + 1) * sum_a;
      double n0 = logspace_add(-c0, 3 * r - c0);
      double n1 = logspace_add(r - c1, 2 * r - c1);
      double lv = n1 - n0;
      if (lv > 27.63) lv = 27.63;
      if (lv < -27.63) lv = -27.63;
      item_logit[j] = lv;
    }
    for (int b = 0; b < static_cast<int>(blocks.size()); ++b) {
      fc_block_logprob(item_logit, blocks[b], block_lp);
      for (int o = 0; o < blocks[b].orr; ++o) {
        table[base + offsets[b] + o] = block_lp[o];
      }
    }
  }
  return table;
}

} // namespace

// [[Rcpp::export]]
List cpp_gibbs_mirt_theta(NumericMatrix theta_, NumericMatrix par_,
                          IntegerMatrix response_, NumericVector theta_mu_,
                          NumericMatrix chol_corr_, double log_diag_sum,
                          double step = 15.0, double lower = -6.0,
                          double upper = 6.0) {
  int N = theta_.nrow();
  int D = theta_.ncol();
  int I = par_.nrow();
  Map<MatrixXd> theta(theta_.begin(), N, D);
  Map<MatrixXd> par(par_.begin(), I, par_.ncol());
  Map<MatrixXd> cholC(chol_corr_.begin(), D, D);
  VectorXd mu(D);
  for (int d = 0; d < D; ++d) mu[d] = theta_mu_[d];

  ThetaGrid grid = build_theta_grid(D, lower, upper, step, mu, cholC,
                                    log_diag_sum);
  BinaryLogProbTable prob_table = mirt_logprob_table(grid, par, D, I);
  std::vector<double> logw(grid.G);
  int sampled = 0;
  for (int person = 0; person < N; ++person) {
    for (int g = 0; g < grid.G; ++g) {
      const size_t base = static_cast<size_t>(g) * I;
      double val = grid.log_prior[g];
      for (int j = 0; j < I; ++j) {
        val += (response_(person, j) == 1) ?
          prob_table.logp1[base + j] : prob_table.logp0[base + j];
      }
      logw[g] = val;
    }
    int idx = sample_log_weight(logw);
    if (idx >= 0) {
      set_theta_row(theta, person, grid, idx);
      ++sampled;
    }
  }
  return List::create(_["theta"] = theta_, _["accepted"] = sampled,
                      _["proposed"] = N);
}

// [[Rcpp::export]]
List cpp_gibbs_mgpcm_theta(NumericMatrix theta_, NumericMatrix par_,
                           IntegerMatrix response_,
                           IntegerVector length_poly_,
                           NumericVector theta_mu_,
                           NumericMatrix chol_corr_, double log_diag_sum,
                           double step = 15.0, double lower = -6.0,
                           double upper = 6.0) {
  int N = theta_.nrow();
  int D = theta_.ncol();
  int I = par_.nrow();
  Map<MatrixXd> theta(theta_.begin(), N, D);
  Map<MatrixXd> par(par_.begin(), I, par_.ncol());
  Map<MatrixXd> cholC(chol_corr_.begin(), D, D);
  VectorXd mu(D);
  for (int d = 0; d < D; ++d) mu[d] = theta_mu_[d];

  ThetaGrid grid = build_theta_grid(D, lower, upper, step, mu, cholC,
                                    log_diag_sum);
  std::vector<int> offsets = category_offsets(length_poly_, I);
  const int total_cat = offsets[I];
  std::vector<double> prob_table = mgpcm_logprob_table(
    grid, par, length_poly_, offsets, total_cat, D, I
  );
  std::vector<double> logw(grid.G);
  int sampled = 0;
  for (int person = 0; person < N; ++person) {
    for (int g = 0; g < grid.G; ++g) {
      const size_t base = static_cast<size_t>(g) * total_cat;
      double val = grid.log_prior[g];
      for (int j = 0; j < I; ++j) {
        int y = response_(person, j);
        int Ki = length_poly_[j];
        if (y >= 0 && y < Ki) val += prob_table[base + offsets[j] + y];
      }
      logw[g] = val;
    }
    int idx = sample_log_weight(logw);
    if (idx >= 0) {
      set_theta_row(theta, person, grid, idx);
      ++sampled;
    }
  }
  return List::create(_["theta"] = theta_, _["accepted"] = sampled,
                      _["proposed"] = N);
}

// [[Rcpp::export]]
List cpp_gibbs_mggum_theta(NumericMatrix theta_, NumericMatrix par_,
                           IntegerMatrix response_,
                           IntegerVector length_poly_,
                           NumericVector theta_mu_,
                           NumericMatrix chol_corr_, double log_diag_sum,
                           double step = 15.0, double lower = -6.0,
                           double upper = 6.0) {
  int N = theta_.nrow();
  int D = theta_.ncol();
  int I = par_.nrow();
  Map<MatrixXd> theta(theta_.begin(), N, D);
  Map<MatrixXd> par(par_.begin(), I, par_.ncol());
  Map<MatrixXd> cholC(chol_corr_.begin(), D, D);
  VectorXd mu(D);
  for (int d = 0; d < D; ++d) mu[d] = theta_mu_[d];

  ThetaGrid grid = build_theta_grid(D, lower, upper, step, mu, cholC,
                                    log_diag_sum);
  std::vector<int> offsets = category_offsets(length_poly_, I);
  const int total_cat = offsets[I];
  std::vector<double> prob_table = mggum_logprob_table(
    grid, par, length_poly_, offsets, total_cat, D, I
  );
  std::vector<double> logw(grid.G);
  int sampled = 0;
  for (int person = 0; person < N; ++person) {
    for (int g = 0; g < grid.G; ++g) {
      const size_t base = static_cast<size_t>(g) * total_cat;
      double val = grid.log_prior[g];
      for (int j = 0; j < I; ++j) {
        int y = response_(person, j);
        int Ki = length_poly_[j];
        if (y >= 0 && y < Ki) val += prob_table[base + offsets[j] + y];
      }
      logw[g] = val;
    }
    int idx = sample_log_weight(logw);
    if (idx >= 0) {
      set_theta_row(theta, person, grid, idx);
      ++sampled;
    }
  }
  return List::create(_["theta"] = theta_, _["accepted"] = sampled,
                      _["proposed"] = N);
}

// [[Rcpp::export]]
List cpp_gibbs_fcmirt_theta(NumericMatrix theta_, NumericMatrix par_,
                            IntegerMatrix response_, NumericVector theta_mu_,
                            NumericMatrix chol_corr_, double log_diag_sum,
                            List patterns_total, List patterns,
                            double step = 15.0, double lower = -6.0,
                            double upper = 6.0) {
  int N = theta_.nrow();
  int D = theta_.ncol();
  int I = par_.nrow();
  int n_blocks = response_.ncol();
  Map<MatrixXd> theta(theta_.begin(), N, D);
  Map<MatrixXd> par(par_.begin(), I, par_.ncol());
  Map<MatrixXd> cholC(chol_corr_.begin(), D, D);
  VectorXd mu(D);
  for (int d = 0; d < D; ++d) mu[d] = theta_mu_[d];

  std::vector<FcBlock> blocks(n_blocks);
  for (int b = 0; b < n_blocks; ++b) {
    blocks[b] = parse_fc(patterns_total[b], patterns[b]);
  }

  ThetaGrid grid = build_theta_grid(D, lower, upper, step, mu, cholC,
                                    log_diag_sum);
  std::vector<int> offsets = fc_response_offsets(blocks);
  const int total_resp = offsets[n_blocks];
  std::vector<double> prob_table = fcmirt_logprob_table(
    grid, par, blocks, offsets, total_resp, D, I
  );
  std::vector<double> logw(grid.G);
  int sampled = 0;
  for (int person = 0; person < N; ++person) {
    for (int g = 0; g < grid.G; ++g) {
      const size_t base = static_cast<size_t>(g) * total_resp;
      double val = grid.log_prior[g];
      for (int b = 0; b < n_blocks; ++b) {
        int y = response_(person, b);
        if (y >= 1 && y <= blocks[b].orr) {
          val += prob_table[base + offsets[b] + y - 1];
        }
      }
      logw[g] = val;
    }
    int idx = sample_log_weight(logw);
    if (idx >= 0) {
      set_theta_row(theta, person, grid, idx);
      ++sampled;
    }
  }
  return List::create(_["theta"] = theta_, _["accepted"] = sampled,
                      _["proposed"] = N);
}

// [[Rcpp::export]]
List cpp_gibbs_fcggum_theta(NumericMatrix theta_, NumericMatrix par_,
                            IntegerMatrix response_, NumericVector theta_mu_,
                            NumericMatrix chol_corr_, double log_diag_sum,
                            List patterns_total, List patterns,
                            double step = 15.0, double lower = -6.0,
                            double upper = 6.0) {
  int N = theta_.nrow();
  int D = theta_.ncol();
  int I = par_.nrow();
  int n_blocks = response_.ncol();
  Map<MatrixXd> theta(theta_.begin(), N, D);
  Map<MatrixXd> par(par_.begin(), I, par_.ncol());
  Map<MatrixXd> cholC(chol_corr_.begin(), D, D);
  VectorXd mu(D);
  for (int d = 0; d < D; ++d) mu[d] = theta_mu_[d];

  std::vector<FcBlock> blocks(n_blocks);
  for (int b = 0; b < n_blocks; ++b) {
    blocks[b] = parse_fc(patterns_total[b], patterns[b]);
  }

  ThetaGrid grid = build_theta_grid(D, lower, upper, step, mu, cholC,
                                    log_diag_sum);
  std::vector<int> offsets = fc_response_offsets(blocks);
  const int total_resp = offsets[n_blocks];
  std::vector<double> prob_table = fcggum_logprob_table(
    grid, par, blocks, offsets, total_resp, D, I
  );
  std::vector<double> logw(grid.G);
  int sampled = 0;
  for (int person = 0; person < N; ++person) {
    for (int g = 0; g < grid.G; ++g) {
      const size_t base = static_cast<size_t>(g) * total_resp;
      double val = grid.log_prior[g];
      for (int b = 0; b < n_blocks; ++b) {
        int y = response_(person, b);
        if (y >= 1 && y <= blocks[b].orr) {
          val += prob_table[base + offsets[b] + y - 1];
        }
      }
      logw[g] = val;
    }
    int idx = sample_log_weight(logw);
    if (idx >= 0) {
      set_theta_row(theta, person, grid, idx);
      ++sampled;
    }
  }
  return List::create(_["theta"] = theta_, _["accepted"] = sampled,
                      _["proposed"] = N);
}
