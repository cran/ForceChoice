// MGGUM Model (distance-based formulation)
// r_i(theta) = sqrt( sum_d a_{id}^2 * (theta_d - delta_{id})^2 )
// psi_{ik} = tau_{ik} * sum_d a_{id}
// C_i = K_i - 1  (number of subjective response categories minus 1)
// M_i = 2 * C_i + 1 = 2 * K_i - 1
//
// Q_matrix: I x D matrix with values -1, 0, 1
//   -1: active dimension with delta <= 0
//    1: active dimension with delta >= 0
//    0: both a and delta fixed to 0 (inactive dimension)
// The Q sign controls the side of delta, not the sign of discrimination.
//
// Identification:
//   - theta_raw ~ std_normal() with theta_mu = 0 (location + scale of theta)
//   - delta side constraints (Q=1: delta >= 0, Q=-1: delta <= 0)
//   - tau ordered constraint via positive threshold gaps:
//       tau_gap > 0 is flattened across items, then accumulated from the
//       highest category downward.
//     Guarantees: tau[i,2] < tau[i,3] < ... < tau[i,K_i] < 0
//
// Default priors (overridable via control):
//   a_free    ~ lognormal(a_mu, a_sigma)              default: lognormal(0, 0.5)
//   delta_pos ~ normal(delta_pos_mu, delta_pos_sigma) default: N(1, 0.5)
//   delta_neg ~ normal(delta_neg_mu, delta_neg_sigma) default: N(-1, 0.5)
//   -tau[,k]  ~ lognormal(tau_mu, tau_sigma) for k >= 2 default: lognormal(0, 0.5)
//   theta     ~ MVN(0, Corr)   Corr ~ lkj_corr_cholesky(1)
data {
  int<lower=1> N_obs;
  int<lower=1> N;
  int<lower=1> I;
  int<lower=1> D;
  int<lower=1> max_poly;

  array[N_obs] int<lower=0, upper=max_poly-1> y;
  array[N_obs] int<lower=1, upper=N> pid;
  array[N_obs] int<lower=1, upper=I> iid;
  array[I] int<lower=2, upper=max_poly> length_poly;

  matrix[I, D] Q_matrix;  // -1/1 = active delta side; 0 = inactive

  // Hyperparameters for a (discrimination): lognormal prior
  real a_mu;               // default 0 -> median=1
  real<lower=0> a_sigma;   // default 0.5 (lognormal scale)

  // Hyperparameters for delta (location): normal prior by sign
  real delta_pos_mu;               // default  1
  real<lower=0> delta_pos_sigma;   // default 0.5
  real delta_neg_mu;               // default -1
  real<lower=0> delta_neg_sigma;   // default 0.5

  // Hyperparameters for tau gaps: lognormal prior
  real tau_mu;             // default 0 -> median=1
  real<lower=0> tau_sigma; // default 0.5 (lognormal scale)

  vector[D] theta_mu;
}

transformed data {
  int a_free_size = 0;
  int delta_pos_size = 0;
  int delta_neg_size = 0;
  array[I] int tau_len;
  array[I + 1] int tau_start;
  int tau_free_total = 0;

  for (i in 1:I) {
    for (j in 1:D) {
      if (Q_matrix[i, j] != 0) {
        a_free_size += 1;
        if (Q_matrix[i, j] == 1) {
          delta_pos_size += 1;
        } else {
          delta_neg_size += 1;
        }
      }
    }
  }

  // Flatten item-specific tau gaps into one vector. Item i uses indices
  // tau_start[i]:(tau_start[i + 1] - 1).
  tau_start[1] = 1;
  for (i in 1:I) {
    tau_len[i] = length_poly[i] - 1;
    tau_free_total += tau_len[i];
    tau_start[i + 1] = tau_start[i] + tau_len[i];
  }
}

parameters {
  vector<lower=0>[a_free_size] a_free;       // positivity constraint
  vector<upper=0>[delta_neg_size] delta_neg; // Q=-1: negative-side delta
  vector<lower=0>[delta_pos_size] delta_pos; // Q=1: positive-side delta
  vector<lower=0>[tau_free_total] tau_gap;   // positive gaps that induce ordered tau values

  cholesky_factor_corr[D] L_Corr;
  array[N] vector[D] theta_raw;
}

transformed parameters {
  matrix<lower=0>[I, D] a; // lower=0 allows zero entries for Q=0
  matrix[I, D] delta;
  matrix[I, max_poly] tau;
  matrix[N, D] theta;

  // Fill a matrix using Q_matrix.
  {
    int idx = 1;
    for (i in 1:I) {
      for (j in 1:D) {
        if (Q_matrix[i, j] != 0) {
          a[i, j] = a_free[idx];
          idx += 1;
        } else {
          a[i, j] = 0.0;
        }
      }
    }
  }

  // Fill delta matrix.
  {
    int idx_pos = 1;
    int idx_neg = 1;
    for (i in 1:I) {
      for (j in 1:D) {
        if (Q_matrix[i, j] == 1) {
          delta[i, j] = delta_pos[idx_pos];
          idx_pos += 1;
        } else if (Q_matrix[i, j] == -1) {
          delta[i, j] = delta_neg[idx_neg];
          idx_neg += 1;
        } else {
          delta[i, j] = 0.0;
        }
      }
    }
  }

  // Fill tau from positive gaps.
  // For Ki = 5 and gaps g1..g4:
  // tau5 = -g4, tau4 = -(g3+g4), ..., tau2 = -(g1+...+g4).
  for (i in 1:I) {
    int Ki = length_poly[i];
    real tau_abs = 0.0;

    tau[i, 1] = 0.0;

    for (k_rev in 1:(Ki - 1)) {
      int k = Ki - k_rev + 1;
      int gap_idx = tau_start[i + 1] - k_rev;
      tau_abs += tau_gap[gap_idx];
      tau[i, k] = -tau_abs;
    }

    if (Ki < max_poly) {
      for (k in (Ki + 1):max_poly) {
        tau[i, k] = 0.0;
      }
    }
  }

  // Compute theta.
  for (n in 1:N) {
    if (D > 1)
      theta[n] = (theta_mu + L_Corr * theta_raw[n])';
    else
      theta[n, 1] = theta_mu[1] + theta_raw[n][1];
  }
}

model {
  vector[I] a_sum;
  matrix[I, max_poly] cum_psi_item;

  // Priors.
  a_free    ~ lognormal(a_mu, a_sigma);
  delta_neg ~ normal(delta_neg_mu, delta_neg_sigma);
  delta_pos ~ normal(delta_pos_mu, delta_pos_sigma);
  for (i in 1:I) {
    for (k in 2:length_poly[i]) {
      real tau_abs = -tau[i, k];
      tau_abs ~ lognormal(tau_mu, tau_sigma);
    }
  }

  if (D > 1)
    L_Corr ~ lkj_corr_cholesky(1);

  for (n in 1:N)
    theta_raw[n] ~ std_normal();

  for (i in 1:I) {
    real cum_psi = 0.0;
    a_sum[i] = sum(a[i]);

    for (k in 1:length_poly[i]) {
      cum_psi += tau[i, k] * a_sum[i];
      cum_psi_item[i, k] = cum_psi;
    }
  }

  // GGUM likelihood (distance-based formulation).
  for (n in 1:N_obs) {
    int i = iid[n];
    int K_i = length_poly[i];
    int M_i = 2 * K_i - 1;

    real r_sq = 0;
    for (d in 1:D)
      r_sq += square(a[i, d]) * square(theta[pid[n], d] - delta[i, d]);
    real r = sqrt(r_sq + 1e-12);

    vector[K_i] log_prob;
    for (k in 1:K_i) {
      log_prob[k] = log_sum_exp(
        (k - 1) * r - cum_psi_item[i, k],
        (M_i - k + 1) * r - cum_psi_item[i, k]
      );
    }

    target += log_prob[y[n] + 1] - log_sum_exp(log_prob);
  }
}

generated quantities {
  matrix[D, D] Corr;
  vector[N_obs] log_lik;

  if (D > 1)
    Corr = multiply_lower_tri_self_transpose(L_Corr);
  else
    Corr[1, 1] = 1.0;

  {
    vector[I] a_sum;
    matrix[I, max_poly] cum_psi_item;

    for (i in 1:I) {
      real cum_psi = 0.0;
      a_sum[i] = sum(a[i]);
      for (k in 1:length_poly[i]) {
        cum_psi += tau[i, k] * a_sum[i];
        cum_psi_item[i, k] = cum_psi;
      }
    }

    for (n in 1:N_obs) {
      int i = iid[n];
      int K_i = length_poly[i];
      int M_i = 2 * K_i - 1;
      real r_sq = 0;
      for (d in 1:D)
        r_sq += square(a[i, d]) * square(theta[pid[n], d] - delta[i, d]);
      real r = sqrt(r_sq + 1e-12);

      vector[K_i] lp;
      for (k in 1:K_i)
        lp[k] = log_sum_exp(
          (k - 1) * r - cum_psi_item[i, k],
          (M_i - k + 1) * r - cum_psi_item[i, k]
        );
      log_lik[n] = lp[y[n] + 1] - log_sum_exp(lp);
    }
  }
}
