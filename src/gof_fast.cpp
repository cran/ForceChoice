#define ARMA_DONT_USE_BLAS
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace {

constexpr double kProbabilityFloor = 1e-12;

inline double probability_floor(const double value, const char* name) {
  if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
    stop("'%s' must contain finite probabilities in [0, 1].", name);
  }
  return std::min(1.0 - kProbabilityFloor,
                  std::max(kProbabilityFloor, value));
}

} // namespace

inline bool grouped_conflict2(const int a, const int b,
                              const IntegerVector& group) {
  return (a != b) && (group[a] == group[b]);
}

inline bool grouped_conflict3(const int a, const int b, const int c,
                              const IntegerVector& group) {
  return grouped_conflict2(a, b, group) ||
         grouped_conflict2(a, c, group) ||
         grouped_conflict2(b, c, group);
}

inline bool grouped_conflict4(const int a, const int b, const int c, const int d,
                              const IntegerVector& group) {
  return grouped_conflict2(a, b, group) ||
         grouped_conflict2(a, c, group) ||
         grouped_conflict2(a, d, group) ||
         grouped_conflict2(b, c, group) ||
         grouped_conflict2(b, d, group) ||
         grouped_conflict2(c, d, group);
}

// Helper: get EIs / EIs2 from category-probability matrix.
// prob: nq x sum(K)  (column-major; columns grouped by item)
// K   : integer vector, length I, number of categories per item
// Returns a list with EIs (nq x I) and EIs2 (nq x I)
// [[Rcpp::export]]
List cpp_compute_EIs(NumericMatrix prob, IntegerVector K) {
  const int nq = prob.nrow();               // number of quadrature points
  const int I  = K.size();                  // number of items

  NumericMatrix EIs (nq, I);
  NumericMatrix EIs2(nq, I);

  int col_offset = 0;
  for (int i = 0; i < I; ++i) {
    const int Ki = K[i];                    // number of categories for this item
    const int col_start = col_offset;       // first column of this item
    const int col_end   = col_start + Ki;   // one past the last column

    // Scores: 0, 1, ..., Ki-1
    NumericVector scores(Ki);
    for (int k = 0; k < Ki; ++k) scores[k] = (double)k;
    NumericVector scores2(Ki);
    for (int k = 0; k < Ki; ++k) scores2[k] = (double)(k * k);

    for (int q = 0; q < nq; ++q) {
      // row-sum of this item's category probabilities
      double row_sum = 0.0;
      for (int c = col_start; c < col_end; ++c) {
        row_sum += prob(q, c);
      }
      // normalise
      double e1  = 0.0;
      double e2  = 0.0;
      for (int c = col_start; c < col_end; ++c) {
        double p_ic = prob(q, c) / row_sum;
        int    kk   = c - col_start;
        e1 += p_ic * scores[kk];
        e2 += p_ic * scores2[kk];
      }
      EIs (q, i) = e1;
      EIs2(q, i) = e2;
    }
    col_offset += Ki;
  }

  return List::create(Named("EIs")  = EIs,
                      Named("EIs2") = EIs2);
}

// Predicted-moment vector (the inner part of moments.predicted.func).
// prob       : nq x (sum(K) if poly else I)
// pi         : length nq, quadrature weights
// pairs_total: 2 x pairs.n  (1-indexed column indices)
// is_poly    : logical
// K          : length I (only used when is_poly)
// Returns a NumericVector of length I + pairs.n
// [[Rcpp::export]]
NumericVector cpp_moments_predicted(NumericMatrix prob, NumericVector pi,
                                    IntegerMatrix pairs_total,
                                    bool is_poly, IntegerVector K) {
  const int nq   = prob.nrow();
  const int I    = is_poly ? K.size() : prob.ncol();
  // --- build EIs ---
  NumericMatrix EIs(nq, I);
  if (is_poly) {
    int col_offset = 0;
    for (int i = 0; i < I; ++i) {
      const int Ki = K[i];
      const int col_start = col_offset;
      const int col_end   = col_start + Ki;
      // scores 0,1,...,Ki-1
      std::vector<double> scores(Ki);
      for (int k = 0; k < Ki; ++k) scores[k] = (double)k;
      for (int q = 0; q < nq; ++q) {
        double row_sum = 0.0;
        for (int c = col_start; c < col_end; ++c) row_sum += prob(q, c);
        double e1 = 0.0;
        for (int c = col_start; c < col_end; ++c) {
          e1 += (prob(q, c) / row_sum) * scores[c - col_start];
        }
        EIs(q, i) = e1;
      }
      col_offset += Ki;
    }
  } else {
    // binary: EIs == prob (already clipped to [1e-10, 1-1e-10])
    for (int i = 0; i < I; ++i) {
      for (int q = 0; q < nq; ++q) {
        EIs(q, i) = prob(q, i);
      }
    }
  }

  const int pairs_n = pairs_total.ncol();
  const int moments_n = I + pairs_n;
  NumericVector res(moments_n);

  // first-order moments
  for (int i = 0; i < I; ++i) {
    double s = 0.0;
    for (int q = 0; q < nq; ++q) {
      s += EIs(q, i) * pi[q];
    }
    res[i] = s;
  }

  // second-order moments (pair products)
  if (pairs_n > 0) {
    for (int k = 0; k < pairs_n; ++k) {
      const int ii = pairs_total(0, k) - 1;  // convert to 0-based
      const int jj = pairs_total(1, k) - 1;
      double s = 0.0;
      for (int q = 0; q < nq; ++q) {
        s += pi[q] * EIs(q, ii) * EIs(q, jj);
      }
      res[I + k] = s;
    }
  }

  return res;
}

// [[Rcpp::export]]
NumericVector cpp_moments_predicted_grouped(NumericMatrix prob, NumericVector pi,
                                            IntegerMatrix pairs_total,
                                            bool is_poly, IntegerVector K,
                                            IntegerVector group) {
  const int nq = prob.nrow();
  const int I = is_poly ? K.size() : prob.ncol();

  NumericMatrix EIs(nq, I);
  if (is_poly) {
    int col_offset = 0;
    for (int i = 0; i < I; ++i) {
      const int Ki = K[i];
      const int col_start = col_offset;
      const int col_end = col_start + Ki;
      std::vector<double> scores(Ki);
      for (int k = 0; k < Ki; ++k) scores[k] = (double)k;
      for (int q = 0; q < nq; ++q) {
        double row_sum = 0.0;
        for (int c = col_start; c < col_end; ++c) row_sum += prob(q, c);
        double e1 = 0.0;
        for (int c = col_start; c < col_end; ++c) {
          e1 += (prob(q, c) / row_sum) * scores[c - col_start];
        }
        EIs(q, i) = e1;
      }
      col_offset += Ki;
    }
  } else {
    for (int i = 0; i < I; ++i) {
      for (int q = 0; q < nq; ++q) {
        EIs(q, i) = prob(q, i);
      }
    }
  }

  const int pairs_n = pairs_total.ncol();
  NumericVector res(I + pairs_n);

  for (int i = 0; i < I; ++i) {
    double s = 0.0;
    for (int q = 0; q < nq; ++q) {
      s += EIs(q, i) * pi[q];
    }
    res[i] = s;
  }

  if (pairs_n > 0) {
    for (int k = 0; k < pairs_n; ++k) {
      const int ii = pairs_total(0, k) - 1;
      const int jj = pairs_total(1, k) - 1;
      double s = 0.0;
      if (!grouped_conflict2(ii, jj, group)) {
        for (int q = 0; q < nq; ++q) {
          s += pi[q] * EIs(q, ii) * EIs(q, jj);
        }
      }
      res[I + k] = s;
    }
  }

  return res;
}

// Xi matrix (Xi11, Xi12, Xi22).
// EIs         : nq x I
// EIs2        : nq x I  (squared expected scores: E[X^2])
// pi          : length nq
// pairs_total : 2 x pairs.n  (1-indexed column indices)
// Returns list(Xi11, Xi12, Xi22, Xi2)
// [[Rcpp::export]]
List cpp_compute_Xi(NumericMatrix EIs, NumericMatrix EIs2, NumericVector pi,
                    IntegerMatrix pairs_total) {
  // --- work in arma for speed ---
  const int nq = EIs.nrow();
  const int I  = EIs.ncol();
  const int pairs_n = pairs_total.ncol();

  // copy to arma (read-only)
  arma::mat  EIs_mat (EIs.begin(),  nq, I,  false);
  arma::mat  EIs2_mat(EIs2.begin(), nq, I,  false);
  arma::vec  pi_vec(pi.begin(), nq, false);

  // first-moment expectations
  arma::vec pa_vec(I);
  for (int i = 0; i < I; ++i) {
    pa_vec(i) = arma::dot(EIs_mat.col(i), pi_vec);
  }

  // === Xi11 ===
  arma::mat Xi11(I, I, arma::fill::zeros);
  for (int i = 0; i < I; ++i) {
    for (int j = 0; j <= i; ++j) {
      double pab;
      if (i == j) {
        pab = arma::dot(EIs2_mat.col(i), pi_vec);
      } else {
        pab = arma::accu(EIs_mat.col(i) % EIs_mat.col(j) % pi_vec);
      }
      double val = pab - pa_vec(i) * pa_vec(j);
      Xi11(i, j) = val;
      Xi11(j, i) = val;
    }
  }

  // === Xi12 ===
  arma::mat Xi12(I, pairs_n, arma::fill::zeros);
  if (pairs_n > 0) {
    for (int k = 0; k < I; ++k) {
      for (int p_idx = 0; p_idx < pairs_n; ++p_idx) {
        const int ii = pairs_total(0, p_idx) - 1;
        const int jj = pairs_total(1, p_idx) - 1;

        double pab  = arma::accu(EIs_mat.col(ii) % EIs_mat.col(jj) % pi_vec);
        double pc   = pa_vec(k);
        double pabc;

        if (ii == k) {
          pabc = arma::accu(EIs2_mat.col(ii) % EIs_mat.col(jj) % pi_vec);
        } else if (jj == k) {
          pabc = arma::accu(EIs_mat.col(ii) % EIs2_mat.col(jj) % pi_vec);
        } else {
          pabc = arma::accu(EIs_mat.col(ii) % EIs_mat.col(jj) % EIs_mat.col(k) % pi_vec);
        }
        Xi12(k, p_idx) = pabc - pab * pc;
      }
    }
  }

  // === Xi22 ===
  arma::mat Xi22(pairs_n, pairs_n, arma::fill::zeros);
  if (pairs_n > 0) {
    for (int kl_idx = 0; kl_idx < pairs_n; ++kl_idx) {
      const int kk = pairs_total(0, kl_idx) - 1;
      const int ll = pairs_total(1, kl_idx) - 1;

      for (int ij_idx = 0; ij_idx < pairs_n; ++ij_idx) {
        const int ii = pairs_total(0, ij_idx) - 1;
        const int jj = pairs_total(1, ij_idx) - 1;

        double pab  = arma::accu(EIs_mat.col(ii) % EIs_mat.col(jj) % pi_vec);
        double pcd  = arma::accu(EIs_mat.col(kk) % EIs_mat.col(ll) % pi_vec);
        double pabcd;

        if ((ii == kk) && (jj == ll)) {
          // perfect overlap of both items
          pabcd = arma::accu(EIs2_mat.col(ii) % EIs2_mat.col(jj) % pi_vec);
        } else if (ii == kk) {
          pabcd = arma::accu(EIs2_mat.col(ii) % EIs_mat.col(jj) % EIs_mat.col(ll) % pi_vec);
        } else if (jj == kk) {
          pabcd = arma::accu(EIs_mat.col(ii) % EIs2_mat.col(jj) % EIs_mat.col(ll) % pi_vec);
        } else if (ii == ll) {
          pabcd = arma::accu(EIs2_mat.col(ii) % EIs_mat.col(jj) % EIs_mat.col(kk) % pi_vec);
        } else if (jj == ll) {
          pabcd = arma::accu(EIs_mat.col(ii) % EIs2_mat.col(jj) % EIs_mat.col(kk) % pi_vec);
        } else {
          pabcd = arma::accu(EIs_mat.col(ii) % EIs_mat.col(jj) % EIs_mat.col(kk) % EIs_mat.col(ll) % pi_vec);
        }
        Xi22(kl_idx, ij_idx) = pabcd - pab * pcd;
      }
    }
  }

  // === assemble full Xi2 and return ===
  NumericMatrix Xi2_out(I + pairs_n, I + pairs_n);
  for (int r = 0; r < I; ++r) {
    for (int c = 0; c < I; ++c) {
      Xi2_out(r, c) = Xi11(r, c);
    }
    for (int c = 0; c < pairs_n; ++c) {
      Xi2_out(r, I + c)          = Xi12(r, c);
      Xi2_out(I + c, r)          = Xi12(r, c);
    }
  }
  for (int r = 0; r < pairs_n; ++r) {
    for (int c = 0; c < pairs_n; ++c) {
      Xi2_out(I + r, I + c) = Xi22(r, c);
    }
  }

  // also return the sub-blocks as R matrices (wrap arma)
  NumericMatrix Xi11_out(I, I);
  NumericMatrix Xi12_out(I, pairs_n);
  NumericMatrix Xi22_out(pairs_n, pairs_n);
  for (int r = 0; r < I; ++r) {
    for (int c = 0; c < I; ++c) {
      Xi11_out(r, c) = Xi11(r, c);
    }
  }
  for (int r = 0; r < I; ++r) {
    for (int c = 0; c < pairs_n; ++c) {
      Xi12_out(r, c) = Xi12(r, c);
    }
  }
  for (int r = 0; r < pairs_n; ++r) {
    for (int c = 0; c < pairs_n; ++c) {
      Xi22_out(r, c) = Xi22(r, c);
    }
  }

  return List::create(
    Named("Xi11") = Xi11_out,
    Named("Xi12") = Xi12_out,
    Named("Xi22") = Xi22_out,
    Named("Xi2")  = Xi2_out
  );
}

// [[Rcpp::export]]
List cpp_compute_Xi_grouped(NumericMatrix EIs, NumericMatrix EIs2, NumericVector pi,
                            IntegerMatrix pairs_total, IntegerVector group) {
  const int nq = EIs.nrow();
  const int I  = EIs.ncol();
  const int pairs_n = pairs_total.ncol();

  arma::mat EIs_mat(EIs.begin(), nq, I, false);
  arma::mat EIs2_mat(EIs2.begin(), nq, I, false);
  arma::vec pi_vec(pi.begin(), nq, false);

  arma::vec pa_vec(I);
  for (int i = 0; i < I; ++i) {
    pa_vec(i) = arma::dot(EIs_mat.col(i), pi_vec);
  }

  arma::mat Xi11(I, I, arma::fill::zeros);
  for (int i = 0; i < I; ++i) {
    for (int j = 0; j <= i; ++j) {
      double pab;
      if (i == j) {
        pab = arma::dot(EIs2_mat.col(i), pi_vec);
      } else if (grouped_conflict2(i, j, group)) {
        pab = 0.0;
      } else {
        pab = arma::accu(EIs_mat.col(i) % EIs_mat.col(j) % pi_vec);
      }
      double val = pab - pa_vec(i) * pa_vec(j);
      Xi11(i, j) = val;
      Xi11(j, i) = val;
    }
  }

  arma::mat Xi12(I, pairs_n, arma::fill::zeros);
  if (pairs_n > 0) {
    for (int k = 0; k < I; ++k) {
      for (int p_idx = 0; p_idx < pairs_n; ++p_idx) {
        const int ii = pairs_total(0, p_idx) - 1;
        const int jj = pairs_total(1, p_idx) - 1;

        double pab = 0.0;
        if (!grouped_conflict2(ii, jj, group)) {
          pab = arma::accu(EIs_mat.col(ii) % EIs_mat.col(jj) % pi_vec);
        }
        double pc = pa_vec(k);
        double pabc;

        if (grouped_conflict3(k, ii, jj, group)) {
          pabc = 0.0;
        } else if (ii == k) {
          pabc = arma::accu(EIs2_mat.col(ii) % EIs_mat.col(jj) % pi_vec);
        } else if (jj == k) {
          pabc = arma::accu(EIs_mat.col(ii) % EIs2_mat.col(jj) % pi_vec);
        } else {
          pabc = arma::accu(EIs_mat.col(ii) % EIs_mat.col(jj) % EIs_mat.col(k) % pi_vec);
        }
        Xi12(k, p_idx) = pabc - pab * pc;
      }
    }
  }

  arma::mat Xi22(pairs_n, pairs_n, arma::fill::zeros);
  if (pairs_n > 0) {
    for (int kl_idx = 0; kl_idx < pairs_n; ++kl_idx) {
      const int kk = pairs_total(0, kl_idx) - 1;
      const int ll = pairs_total(1, kl_idx) - 1;

      for (int ij_idx = 0; ij_idx < pairs_n; ++ij_idx) {
        const int ii = pairs_total(0, ij_idx) - 1;
        const int jj = pairs_total(1, ij_idx) - 1;

        double pab = 0.0;
        if (!grouped_conflict2(ii, jj, group)) {
          pab = arma::accu(EIs_mat.col(ii) % EIs_mat.col(jj) % pi_vec);
        }
        double pcd = 0.0;
        if (!grouped_conflict2(kk, ll, group)) {
          pcd = arma::accu(EIs_mat.col(kk) % EIs_mat.col(ll) % pi_vec);
        }

        double pabcd;
        if (grouped_conflict4(ii, jj, kk, ll, group)) {
          pabcd = 0.0;
        } else if ((ii == kk) && (jj == ll)) {
          pabcd = arma::accu(EIs2_mat.col(ii) % EIs2_mat.col(jj) % pi_vec);
        } else if (ii == kk) {
          pabcd = arma::accu(EIs2_mat.col(ii) % EIs_mat.col(jj) % EIs_mat.col(ll) % pi_vec);
        } else if (jj == kk) {
          pabcd = arma::accu(EIs_mat.col(ii) % EIs2_mat.col(jj) % EIs_mat.col(ll) % pi_vec);
        } else if (ii == ll) {
          pabcd = arma::accu(EIs2_mat.col(ii) % EIs_mat.col(jj) % EIs_mat.col(kk) % pi_vec);
        } else if (jj == ll) {
          pabcd = arma::accu(EIs_mat.col(ii) % EIs2_mat.col(jj) % EIs_mat.col(kk) % pi_vec);
        } else {
          pabcd = arma::accu(EIs_mat.col(ii) % EIs_mat.col(jj) % EIs_mat.col(kk) % EIs_mat.col(ll) % pi_vec);
        }
        Xi22(kl_idx, ij_idx) = pabcd - pab * pcd;
      }
    }
  }

  NumericMatrix Xi2_out(I + pairs_n, I + pairs_n);
  for (int r = 0; r < I; ++r) {
    for (int c = 0; c < I; ++c) {
      Xi2_out(r, c) = Xi11(r, c);
    }
    for (int c = 0; c < pairs_n; ++c) {
      Xi2_out(r, I + c) = Xi12(r, c);
      Xi2_out(I + c, r) = Xi12(r, c);
    }
  }
  for (int r = 0; r < pairs_n; ++r) {
    for (int c = 0; c < pairs_n; ++c) {
      Xi2_out(I + r, I + c) = Xi22(r, c);
    }
  }

  NumericMatrix Xi11_out(I, I);
  NumericMatrix Xi12_out(I, pairs_n);
  NumericMatrix Xi22_out(pairs_n, pairs_n);
  for (int r = 0; r < I; ++r) {
    for (int c = 0; c < I; ++c) {
      Xi11_out(r, c) = Xi11(r, c);
    }
  }
  for (int r = 0; r < I; ++r) {
    for (int c = 0; c < pairs_n; ++c) {
      Xi12_out(r, c) = Xi12(r, c);
    }
  }
  for (int r = 0; r < pairs_n; ++r) {
    for (int c = 0; c < pairs_n; ++c) {
      Xi22_out(r, c) = Xi22(r, c);
    }
  }

  return List::create(
    Named("Xi11") = Xi11_out,
    Named("Xi12") = Xi12_out,
    Named("Xi22") = Xi22_out,
    Named("Xi2")  = Xi2_out
  );
}

// E2 matrix (for SRMSR).
// [[Rcpp::export]]
NumericMatrix cpp_compute_E2(NumericMatrix EIs, NumericMatrix EIs2,
                             NumericVector pi) {
  const int nq = EIs.nrow();
  const int I  = EIs.ncol();

  arma::mat  EIs_mat (EIs.begin(),  nq, I, false);
  arma::mat  EIs2_mat(EIs2.begin(), nq, I, false);
  arma::vec  pi_vec(pi.begin(), nq, false);

  // E[X_i^2] expectation column
  arma::vec E11(I);
  for (int i = 0; i < I; ++i) {
    E11(i) = arma::dot(EIs2_mat.col(i), pi_vec);
  }

  NumericMatrix E2_out(I, I);
  for (int i = 0; i < I; ++i) {
    for (int j = 0; j <= i; ++j) {
      double val;
      if (i == j) {
        val = E11(i);
      } else {
        val = arma::accu(EIs_mat.col(i) % EIs_mat.col(j) % pi_vec);
      }
      E2_out(i, j) = val;
      E2_out(j, i) = val;
    }
  }
  return E2_out;
}

// [[Rcpp::export]]
NumericMatrix cpp_compute_E2_grouped(NumericMatrix EIs, NumericMatrix EIs2,
                                     NumericVector pi, IntegerVector group) {
  const int nq = EIs.nrow();
  const int I = EIs.ncol();

  arma::mat EIs_mat(EIs.begin(), nq, I, false);
  arma::mat EIs2_mat(EIs2.begin(), nq, I, false);
  arma::vec pi_vec(pi.begin(), nq, false);

  arma::vec E11(I);
  for (int i = 0; i < I; ++i) {
    E11(i) = arma::dot(EIs2_mat.col(i), pi_vec);
  }

  NumericMatrix E2_out(I, I);
  for (int i = 0; i < I; ++i) {
    for (int j = 0; j <= i; ++j) {
      double val;
      if (i == j) {
        val = E11(i);
      } else if (grouped_conflict2(i, j, group)) {
        val = 0.0;
      } else {
        val = arma::accu(EIs_mat.col(i) % EIs_mat.col(j) % pi_vec);
      }
      E2_out(i, j) = val;
      E2_out(j, i) = val;
    }
  }
  return E2_out;
}

// [[Rcpp::export]]
IntegerMatrix cpp_make_pairs(int n) {
  if (n < 1) {
    stop("'n' must be at least 1.");
  }
  const long long n_pairs_ll = static_cast<long long>(n) * (n - 1) / 2;
  if (n_pairs_ll > std::numeric_limits<int>::max()) {
    stop("The number of pairs is too large to allocate.");
  }
  const int n_pairs = static_cast<int>(n_pairs_ll);
  IntegerMatrix out(2, n_pairs);
  int idx = 0;
  for (int i = 1; i < n; ++i) {
    for (int j = i + 1; j <= n; ++j) {
      out(0, idx) = i;
      out(1, idx) = j;
      ++idx;
    }
  }
  return out;
}

// [[Rcpp::export]]
IntegerMatrix cpp_filter_pairs_between_groups(IntegerMatrix pairs_total,
                                              IntegerVector group) {
  if (pairs_total.nrow() != 2) {
    stop("'pairs_total' must have two rows.");
  }

  std::vector<int> keep;
  keep.reserve(pairs_total.ncol());
  for (int k = 0; k < pairs_total.ncol(); ++k) {
    const int ii = pairs_total(0, k);
    const int jj = pairs_total(1, k);
    if (ii < 1 || jj < 1 || ii > group.size() || jj > group.size()) {
      stop("'pairs_total' contains indices outside 'group'.");
    }
    if (group[ii - 1] != group[jj - 1]) {
      keep.push_back(k);
    }
  }

  IntegerMatrix out(2, static_cast<int>(keep.size()));
  for (int k = 0; k < static_cast<int>(keep.size()); ++k) {
    out(0, k) = pairs_total(0, keep[k]);
    out(1, k) = pairs_total(1, keep[k]);
  }
  return out;
}

// [[Rcpp::export]]
List cpp_gof_observed_summary(NumericMatrix response,
                              IntegerMatrix pairs_total) {
  const int N = response.nrow();
  const int I = response.ncol();
  if (pairs_total.nrow() != 2) {
    stop("'pairs_total' must have two rows.");
  }

  const int pairs_n = pairs_total.ncol();
  NumericVector means(I);
  NumericVector seconds(I);
  NumericVector moments(I + pairs_n);
  NumericMatrix cross(I, I);
  NumericMatrix count(I, I);

  for (int i = 0; i < I; ++i) {
    double sum_i = 0.0;
    double sumsq_i = 0.0;
    int n_i = 0;
    for (int r = 0; r < N; ++r) {
      const double xi = response(r, i);
      if (!NumericMatrix::is_na(xi)) {
        sum_i += xi;
        sumsq_i += xi * xi;
        ++n_i;
      }
    }
    means[i] = n_i > 0 ? sum_i / static_cast<double>(n_i) : R_NaN;
    seconds[i] = n_i > 0 ? sumsq_i / static_cast<double>(n_i) : R_NaN;
    moments[i] = means[i];
  }

  for (int i = 0; i < I; ++i) {
    for (int j = i; j < I; ++j) {
      double s = 0.0;
      int n_ij = 0;
      for (int r = 0; r < N; ++r) {
        const double xi = response(r, i);
        const double xj = response(r, j);
        if (!NumericMatrix::is_na(xi) && !NumericMatrix::is_na(xj)) {
          s += xi * xj;
          ++n_ij;
        }
      }
      cross(i, j) = s;
      cross(j, i) = s;
      count(i, j) = n_ij;
      count(j, i) = n_ij;
    }
  }

  for (int k = 0; k < pairs_n; ++k) {
    const int ii = pairs_total(0, k) - 1;
    const int jj = pairs_total(1, k) - 1;
    if (ii < 0 || jj < 0 || ii >= I || jj >= I) {
      stop("'pairs_total' contains invalid response-column indices.");
    }
    moments[I + k] = count(ii, jj) > 0.0 ? cross(ii, jj) / count(ii, jj) : R_NaN;
  }

  return List::create(
    Named("moments") = moments,
    Named("means") = means,
    Named("seconds") = seconds,
    Named("cross") = cross,
    Named("count") = count
  );
}

// [[Rcpp::export]]
double cpp_gof_null_loglik(NumericMatrix response,
                           IntegerVector K,
                           std::string response_type,
                           IntegerVector nominal_groups,
                           NumericVector first_moments) {
  const int N = response.nrow();
  const int I = response.ncol();
  double out = 0.0;

  if (response_type == "polytomous") {
    if (K.size() != I) {
      stop("'K' must have one entry per response column.");
    }
    for (int i = 0; i < I; ++i) {
      const int Ki = K[i];
      std::vector<int> tab(Ki, 0);
      int n_i = 0;
      for (int r = 0; r < N; ++r) {
        const double yi = response(r, i);
        if (NumericMatrix::is_na(yi)) {
          continue;
        }
        const int cat = static_cast<int>(yi);
        if (cat >= 0 && cat < Ki) {
          ++tab[cat];
          ++n_i;
        }
      }
      if (n_i == 0) {
        continue;
      }
      for (int k = 0; k < Ki; ++k) {
        if (tab[k] > 0) {
          out += static_cast<double>(tab[k]) *
            std::log(std::max(kProbabilityFloor,
                              static_cast<double>(tab[k]) / static_cast<double>(n_i)));
        }
      }
    }
    return out;
  }

  if (response_type == "nominal") {
    if (nominal_groups.size() != I) {
      stop("'nominal_groups' must have one entry per response column.");
    }

    std::vector<int> groups;
    groups.reserve(I);
    for (int i = 0; i < I; ++i) {
      if (std::find(groups.begin(), groups.end(), nominal_groups[i]) == groups.end()) {
        groups.push_back(nominal_groups[i]);
      }
    }

    for (std::size_t gidx = 0; gidx < groups.size(); ++gidx) {
      std::vector<int> cols;
      for (int i = 0; i < I; ++i) {
        if (nominal_groups[i] == groups[gidx]) {
          cols.push_back(i);
        }
      }

      std::vector<double> p_cat(cols.size(), kProbabilityFloor);
      double p_sum = 0.0;
      for (std::size_t c = 0; c < cols.size(); ++c) {
        double s = 0.0;
        int n_c = 0;
        for (int r = 0; r < N; ++r) {
          const double y = response(r, cols[c]);
          if (!NumericMatrix::is_na(y)) {
            s += y;
            ++n_c;
          }
        }
        p_cat[c] = std::max(kProbabilityFloor,
                            n_c > 0 ? s / static_cast<double>(n_c) : 0.0);
        p_sum += p_cat[c];
      }
      const double p_ref = std::max(kProbabilityFloor, 1.0 - p_sum);

      for (int r = 0; r < N; ++r) {
        bool all_na = true;
        int hit = -1;
        for (std::size_t c = 0; c < cols.size(); ++c) {
          const double y = response(r, cols[c]);
          if (!NumericMatrix::is_na(y)) {
            all_na = false;
            if (hit < 0 && y == 1.0) {
              hit = static_cast<int>(c);
            }
          }
        }
        if (all_na) {
          continue;
        }
        out += hit >= 0 ? std::log(p_cat[hit]) : std::log(p_ref);
      }
    }
    return out;
  }

  if (first_moments.size() != I) {
    stop("'first_moments' must have one entry per response column.");
  }
  for (int i = 0; i < I; ++i) {
    const double pr = probability_floor(first_moments[i], "first_moments");
    const double log_pr = std::log(pr);
    const double log_qr = std::log1p(-pr);
    for (int r = 0; r < N; ++r) {
      const double y = response(r, i);
      if (!NumericMatrix::is_na(y)) {
        if (y == 1.0) {
          out += log_pr;
        } else if (y == 0.0) {
          out += log_qr;
        } else {
          out += y * log_pr + (1.0 - y) * log_qr;
        }
      }
    }
  }
  return out;
}
