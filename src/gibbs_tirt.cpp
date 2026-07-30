// TIRT block Gibbs theta sampler.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "istem_gibbs_grid.h"
#include <cmath>
#include <vector>

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

using forcechoice_istem::ThetaGrid;
using forcechoice_istem::build_theta_grid;
using forcechoice_istem::sample_log_weight;

template <typename MatrixType>
void set_theta_row(MatrixType& theta, int person, const ThetaGrid& grid,
                   int grid_index) {
  const double* x = &grid.point[static_cast<size_t>(grid_index) * grid.D];
  for (int d = 0; d < grid.D; ++d) theta(person, d) = x[d];
}

double tirt_person_loglik(const double* x,
                          const MatrixXd& Q,
                          const VectorXd& lambda,
                          const VectorXd& psi,
                          const VectorXd& gamma,
                          const std::vector<int>& obs,
                          const std::vector<int>& item_i,
                          const std::vector<int>& item_k,
                          const std::vector<int>& pair_idx,
                          const std::vector<int>& y,
                          int D,
                          int n_items) {
  std::vector<double> qtheta(n_items, 0.0);
  for (int item = 0; item < n_items; ++item) {
    double val = 0.0;
    for (int d = 0; d < D; ++d) val += Q(item, d) * x[d];
    qtheta[item] = val;
  }

  double out = 0.0;
  for (int obs_idx : obs) {
    int i = item_i[obs_idx];
    int k = item_k[obs_idx];
    int pair = pair_idx[obs_idx];
    double mean_diff = lambda[i] * qtheta[i] - lambda[k] * qtheta[k] -
      gamma[pair];
    double denom = std::sqrt(psi[i] * psi[i] + psi[k] * psi[k]);
    double p = R::pnorm(mean_diff / denom, 0.0, 1.0, 1, 0);
    if (p < 1e-12) p = 1e-12;
    if (p > 1.0 - 1e-12) p = 1.0 - 1e-12;
    out += (y[obs_idx] == 1) ? std::log(p) : std::log1p(-p);
  }
  return out;
}

std::vector<double> tirt_logprob_table(const ThetaGrid& grid,
                                       const MatrixXd& Q,
                                       const VectorXd& lambda,
                                       const VectorXd& psi,
                                       const VectorXd& gamma,
                                       const std::vector<int>& item_i,
                                       const std::vector<int>& item_k,
                                       const std::vector<int>& pair_idx,
                                       const std::vector<int>& y,
                                       int D,
                                       int n_items,
                                       int n_obs) {
  std::vector<double> table(static_cast<size_t>(grid.G) * n_obs);
  std::vector<double> qtheta(n_items, 0.0);
  for (int g = 0; g < grid.G; ++g) {
    const double* x = &grid.point[static_cast<size_t>(g) * D];
    const size_t base = static_cast<size_t>(g) * n_obs;
    for (int item = 0; item < n_items; ++item) {
      double val = 0.0;
      for (int d = 0; d < D; ++d) val += Q(item, d) * x[d];
      qtheta[item] = val;
    }

    for (int obs = 0; obs < n_obs; ++obs) {
      int i = item_i[obs];
      int k = item_k[obs];
      int pair = pair_idx[obs];
      double mean_diff = lambda[i] * qtheta[i] - lambda[k] * qtheta[k] -
        gamma[pair];
      double denom = std::sqrt(psi[i] * psi[i] + psi[k] * psi[k]);
      double p = R::pnorm(mean_diff / denom, 0.0, 1.0, 1, 0);
      if (p < 1e-12) p = 1e-12;
      if (p > 1.0 - 1e-12) p = 1.0 - 1e-12;
      table[base + obs] = (y[obs] == 1) ? std::log(p) : std::log1p(-p);
    }
  }
  return table;
}

} // namespace

// [[Rcpp::export]]
List cpp_gibbs_tirt_theta(
    NumericMatrix theta_, NumericMatrix Q_matrix_,
    NumericVector lambda_, NumericVector psi_, NumericVector gamma_,
    IntegerVector Y_, List obs_by_person_,
    IntegerVector Idx_item_i_, IntegerVector Idx_item_k_,
    IntegerVector Idx_pair_,
    NumericVector theta_mu_, NumericMatrix chol_corr_, double log_diag_sum,
    double step = 15.0, double lower = -6.0, double upper = 6.0
) {
  int N = theta_.nrow();
  int D = theta_.ncol();
  int n_items = Q_matrix_.nrow();
  int n_obs = Y_.size();

  Map<MatrixXd> theta(theta_.begin(), N, D);
  Map<MatrixXd> Q(Q_matrix_.begin(), n_items, D);
  Map<MatrixXd> cholC(chol_corr_.begin(), D, D);

  VectorXd lambda(lambda_.size());
  VectorXd psi(psi_.size());
  VectorXd gamma(gamma_.size());
  VectorXd mu(D);
  for (int i = 0; i < lambda_.size(); ++i) lambda[i] = lambda_[i];
  for (int i = 0; i < psi_.size(); ++i) psi[i] = psi_[i];
  for (int i = 0; i < gamma_.size(); ++i) gamma[i] = gamma_[i];
  for (int d = 0; d < D; ++d) mu[d] = theta_mu_[d];

  std::vector<std::vector<int>> obs_by_person(N);
  for (int p = 0; p < N; ++p) {
    IntegerVector obs = as<IntegerVector>(obs_by_person_[p]);
    obs_by_person[p].resize(obs.size());
    for (int j = 0; j < obs.size(); ++j) obs_by_person[p][j] = obs[j] - 1;
  }

  std::vector<int> item_i(n_obs), item_k(n_obs), pair_idx(n_obs), y(n_obs);
  for (int obs = 0; obs < n_obs; ++obs) {
    item_i[obs] = Idx_item_i_[obs] - 1;
    item_k[obs] = Idx_item_k_[obs] - 1;
    pair_idx[obs] = Idx_pair_[obs] - 1;
    y[obs] = Y_[obs];
  }

  ThetaGrid grid = build_theta_grid(D, lower, upper, step, mu, cholC,
                                    log_diag_sum);
  std::vector<double> prob_table = tirt_logprob_table(
    grid, Q, lambda, psi, gamma, item_i, item_k, pair_idx, y, D, n_items,
    n_obs
  );
  std::vector<double> logw(grid.G);
  int sampled = 0;
  for (int person = 0; person < N; ++person) {
    for (int g = 0; g < grid.G; ++g) {
      const size_t base = static_cast<size_t>(g) * n_obs;
      double val = grid.log_prior[g];
      for (int obs_idx : obs_by_person[person]) {
        if (obs_idx >= 0 && obs_idx < n_obs) val += prob_table[base + obs_idx];
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
