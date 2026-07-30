#ifndef FORCECHOICE_ISTEM_GIBBS_GRID_H
#define FORCECHOICE_ISTEM_GIBBS_GRID_H

#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include <vector>

namespace forcechoice_istem {

constexpr int kMaxGridPoints = 5000000;

inline double logspace_add(double a, double b) {
  if (!std::isfinite(a)) return b;
  if (!std::isfinite(b)) return a;
  if (a > b) return a + std::log1p(std::exp(b - a));
  return b + std::log1p(std::exp(a - b));
}

inline double logspace_sum(const std::vector<double>& x) {
  double out = -std::numeric_limits<double>::infinity();
  for (double v : x) out = logspace_add(out, v);
  return out;
}

inline int grid_length_from_step(double step) {
  if (!std::isfinite(step)) {
    Rcpp::stop("'L' must be finite.");
  }
  int L = static_cast<int>(std::floor(step + 0.5));
  if (L < 2) {
    Rcpp::stop("'L' must be at least 2.");
  }
  return L;
}

inline std::vector<double> make_axis(double lower, double upper, int L) {
  if (!std::isfinite(lower) || !std::isfinite(upper) || lower >= upper) {
    Rcpp::stop("'theta.lower' must be smaller than 'theta.upper'.");
  }
  std::vector<double> axis(L);
  double h = (upper - lower) / static_cast<double>(L - 1);
  for (int i = 0; i < L; ++i) axis[i] = lower + h * i;
  return axis;
}

inline double log_mvn_prior(const double* x, int D, const Eigen::VectorXd& mu,
                            const Eigen::MatrixXd& cholC,
                            double log_diag_sum) {
  std::vector<double> z(D);
  double quad = 0.0;
  for (int i = 0; i < D; ++i) {
    double acc = x[i] - mu[i];
    for (int j = 0; j < i; ++j) acc -= cholC(j, i) * z[j];
    z[i] = acc / cholC(i, i);
    quad += z[i] * z[i];
  }
  return -0.5 * quad - log_diag_sum;
}

struct ThetaGrid {
  int D;
  int L;
  int G;
  std::vector<double> point;
  std::vector<double> log_prior;
};

inline ThetaGrid build_theta_grid(int D, double lower, double upper,
                                  double step, const Eigen::VectorXd& mu,
                                  const Eigen::MatrixXd& cholC,
                                  double log_diag_sum) {
  int L = grid_length_from_step(step);
  long long Gll = 1;
  for (int d = 0; d < D; ++d) {
    Gll *= L;
    if (Gll > kMaxGridPoints) {
      Rcpp::stop("'L'^D is too large for block Gibbs sampling.");
    }
  }

  ThetaGrid grid;
  grid.D = D;
  grid.L = L;
  grid.G = static_cast<int>(Gll);
  grid.point.resize(static_cast<size_t>(grid.G) * D);
  grid.log_prior.resize(grid.G);

  std::vector<double> axis = make_axis(lower, upper, L);
  for (int g = 0; g < grid.G; ++g) {
    int idx = g;
    for (int d = 0; d < D; ++d) {
      grid.point[static_cast<size_t>(g) * D + d] = axis[idx % L];
      idx /= L;
    }
    grid.log_prior[g] = log_mvn_prior(
      &grid.point[static_cast<size_t>(g) * D], D, mu, cholC, log_diag_sum
    );
  }

  return grid;
}

inline int sample_log_weight(const std::vector<double>& logw) {
  double max_log = -std::numeric_limits<double>::infinity();
  for (double v : logw) {
    if (std::isfinite(v) && v > max_log) max_log = v;
  }
  if (!std::isfinite(max_log)) return -1;

  double total = 0.0;
  for (double v : logw) {
    if (std::isfinite(v)) total += std::exp(v - max_log);
  }
  if (!std::isfinite(total) || total <= 0.0) return -1;

  double u = R::runif(0.0, total);
  double acc = 0.0;
  for (int i = 0; i < static_cast<int>(logw.size()); ++i) {
    if (std::isfinite(logw[i])) {
      acc += std::exp(logw[i] - max_log);
      if (u <= acc) return i;
    }
  }
  return static_cast<int>(logw.size()) - 1;
}

} // namespace forcechoice_istem

#endif
