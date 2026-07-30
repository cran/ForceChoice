// Item log-likelihood functions used by R's optim() for MAP item updates.
// These are kept separate from the Gibbs theta samplers.
#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {
constexpr double kMirtFloor = 1e-10;

inline double inv_logit(double x) {
  if (x >= 0.0) { double z = std::exp(-x); return 1.0 / (1.0 + z); }
  double z = std::exp(x); return z / (1.0 + z);
}
inline double logspace_add(double a, double b) {
  if (!std::isfinite(a)) return b; if (!std::isfinite(b)) return a;
  if (a > b) return a + std::log1p(std::exp(b - a));
  return b + std::log1p(std::exp(a - b));
}
double logspace_sum(const std::vector<double>& x) {
  double o = -1e300; for (double v : x) o = logspace_add(o, v); return o;
}
} // namespace

// [[Rcpp::export]]
double cpp_mirt_item_loglik_weighted(NumericMatrix theta, NumericVector a,
                                     double b, double c, double upper,
                                     IntegerVector response,
                                     NumericVector weight) {
  int N = theta.nrow(), D = theta.ncol();
  double out = 0.0;
  if (response.size() != N || weight.size() != N) {
    stop("'response' and 'weight' must have one element per row of theta.");
  }
  for (int p = 0; p < N; ++p) {
    double eta = -b;
    for (int d = 0; d < D; ++d) eta += theta(p, d) * a[d];
    double pr = c + (upper - c) * inv_logit(eta);
    if (pr < kMirtFloor) pr = kMirtFloor;
    if (pr > 1.0 - kMirtFloor) pr = 1.0 - kMirtFloor;
    out += weight[p] * ((response[p] == 1) ? std::log(pr) : std::log1p(-pr));
  }
  return out;
}

// [[Rcpp::export]]
double cpp_mgpcm_item_loglik_weighted(NumericMatrix theta, NumericVector a,
                                      NumericVector d, IntegerVector response,
                                      NumericVector weight) {
  int N = theta.nrow(), D = theta.ncol();
  int Ki = d.size();
  double out = 0.0;
  std::vector<double> theta_i(D);
  if (response.size() != N || weight.size() != N) {
    stop("'response' and 'weight' must have one element per row of theta.");
  }
  for (int p = 0; p < N; ++p) {
    for (int dim = 0; dim < D; ++dim) theta_i[dim] = theta(p, dim);
    double eta = 0.0;
    for (int dim = 0; dim < D; ++dim) eta += theta_i[dim] * a[dim];
    std::vector<double> log_num(Ki);
    for (int k = 0; k < Ki; ++k)
      log_num[k] = static_cast<double>(k) * eta + d[k];
    double ls = logspace_sum(log_num);
    int y = response[p];
    if (y >= 0 && y < Ki) out += weight[p] * (log_num[y] - ls);
  }
  return out;
}

// [[Rcpp::export]]
double cpp_mggum_item_loglik_weighted(NumericMatrix theta, NumericVector par_i,
                                      int Ki, IntegerVector response,
                                      NumericVector weight) {
  int N = theta.nrow(), D = theta.ncol();
  double out = 0.0;
  std::vector<double> theta_i(D);
  if (response.size() != N || weight.size() != N) {
    stop("'response' and 'weight' must have one element per row of theta.");
  }
  for (int p = 0; p < N; ++p) {
    for (int dim = 0; dim < D; ++dim) theta_i[dim] = theta(p, dim);
    double sum_a = 0.0, r_sq = 0.0;
    for (int d = 0; d < D; ++d) {
      sum_a += par_i[d];
      double diff = theta_i[d] - par_i[D + d];
      r_sq += par_i[d] * par_i[d] * diff * diff;
    }
    double r_i = std::sqrt(r_sq + 1e-10);
    int Mi = 2 * Ki - 1;
    std::vector<double> log_num(Ki);
    double cum_psi = 0.0;
    for (int k = 0; k < Ki; ++k) {
      cum_psi += par_i[2 * D + k] * sum_a;
      double z = static_cast<double>(k);
      log_num[k] = logspace_add(z * r_i - cum_psi,
                                (static_cast<double>(Mi) - z) * r_i - cum_psi);
    }
    double ls = logspace_sum(log_num);
    int y = response[p];
    if (y >= 0 && y < Ki) out += weight[p] * (log_num[y] - ls);
  }
  return out;
}

// ============================================================================
// TIRT per-item log-likelihood (for fast per-item optimisation of lambda, psi)
// ============================================================================
// All inputs except lambda_val and psi_val are pre-computed once in R and
// reused across optim() iterations; C++ only evaluates the probit loop.
//
// qt_target[p]   = qtheta[person, target_item]
// qt_other[p]    = qtheta[person, other_item]
// lambda_other   = lambda[other_item]
// psi_other      = psi[other_item]
// gamma_pair     = gamma[pair_index]
// target_is_first = 1 when target == item_i, 0 when target == item_k
// ============================================================================
// [[Rcpp::export]]
double cpp_tirt_item_loglik_weighted(
    NumericVector qt_target,
    NumericVector qt_other,
    NumericVector lambda_other,
    NumericVector psi_other,
    NumericVector gamma_pair,
    IntegerVector target_is_first,
    IntegerVector Y,
    NumericVector weight,
    double lambda_val,
    double psi_val
) {
  int N = Y.size();
  double ll = 0.0;
  for (int p = 0; p < N; ++p) {
    double mu_val, den;
    if (target_is_first[p]) {
      mu_val = lambda_val * qt_target[p]
             - lambda_other[p] * qt_other[p]
             - gamma_pair[p];
    } else {
      mu_val = lambda_other[p] * qt_other[p]
             - lambda_val * qt_target[p]
             - gamma_pair[p];
    }
    den = std::sqrt(psi_val * psi_val + psi_other[p] * psi_other[p]);
    double pr = R::pnorm(mu_val / den, 0.0, 1.0, 1, 0);
    if (pr < 1e-12) pr = 1e-12;
    if (pr > 1.0 - 1e-12) pr = 1.0 - 1e-12;
    ll += weight[p] * ((Y[p] == 1) ? std::log(pr) : std::log1p(-pr));
  }
  return ll;
}

// ============================================================================
// TIRT per-pair log-likelihood (for fast per-pair optimisation of gamma)
// ============================================================================
// mu_no_gamma[p] = lambda[ii]*qtheta[p,ii] - lambda[ik]*qtheta[p,ik]
// denom[p]       = sqrt(psi[ii]^2 + psi[ik]^2)
// Both are pre-computed once in R and reused across optim() iterations.
// ============================================================================
// [[Rcpp::export]]
double cpp_tirt_pair_loglik_weighted(
    NumericVector mu_no_gamma,
    NumericVector denom,
    IntegerVector Y,
    NumericVector weight,
    double gamma_val
) {
  int N = Y.size();
  double ll = 0.0;
  for (int p = 0; p < N; ++p) {
    double pr = R::pnorm((mu_no_gamma[p] - gamma_val) / denom[p], 0.0, 1.0, 1, 0);
    if (pr < 1e-12) pr = 1e-12;
    if (pr > 1.0 - 1e-12) pr = 1.0 - 1e-12;
    ll += weight[p] * ((Y[p] == 1) ? std::log(pr) : std::log1p(-pr));
  }
  return ll;
}
