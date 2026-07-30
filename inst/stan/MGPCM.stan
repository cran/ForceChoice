// MGPCM Model
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

  matrix[I, D] Q_matrix;

  real a_mu;
  real<lower=0> a_sigma;
  real d_mu;
  real<lower=0> d_sigma;

  vector[D] theta_mu;
}

transformed data {
  int a_free_size = 0;
  for (i in 1:I) {
    for (j in 1:D) {
      if (Q_matrix[i, j] == 1) {
        a_free_size += 1;
      }
    }
  }
  int d_free_size = sum(length_poly) - I;
}

parameters {
  vector<lower=0>[a_free_size] a_free;
  vector[d_free_size] d_free;
  cholesky_factor_corr[D] L_Corr;
  array[N] vector[D] theta_raw;
}

transformed parameters {
  matrix<lower=0>[I, D] a;
  matrix[I, max_poly] d;
  matrix[N, D] theta;
  vector[N_obs] eta;

  // Fill a matrix using Q_matrix
  {
    int idx = 1;
    for (i in 1:I) {
      for (j in 1:D) {
        if (Q_matrix[i, j] == 1) {
          a[i, j] = a_free[idx];
          idx += 1;
        } else {
          a[i, j] = 0.0;
        }
      }
    }
  }

  // Fill d matrix: d[i, 1] = 0 (fixed), d[i, 2..length_poly[i]] free
  {
    int idx = 1;
    for (i in 1:I) {
      d[i, 1] = 0.0;
      for (k in 2:length_poly[i]) {
        d[i, k] = d_free[idx];
        idx += 1;
      }
      if (length_poly[i] < max_poly) {
        for (k in (length_poly[i] + 1):max_poly) {
          d[i, k] = 0.0;
        }
      }
    }
  }

  // Compute theta
  for (n in 1:N) {
    if (D > 1) {
      theta[n] = (theta_mu + L_Corr * theta_raw[n])';
    } else {
      theta[n, 1] = theta_mu[1] + theta_raw[n][1];
    }
  }

  // Linear predictor for each observation
  for (n in 1:N_obs) {
    eta[n] = dot_product(theta[pid[n]], a[iid[n]]);
  }
}

model {
  a_free ~ lognormal(a_mu, a_sigma);
  d_free ~ normal(d_mu, d_sigma);

  if (D > 1) {
    L_Corr ~ lkj_corr_cholesky(1);
  }

  for (n in 1:N) {
    theta_raw[n] ~ std_normal();
  }

  // Categorical likelihood (GPCM)
  for (n in 1:N_obs) {
    int i = iid[n];
    int K = length_poly[i];
    vector[K] log_prob;

    for (k in 1:K) {
      log_prob[k] = (k - 1) * eta[n] + d[i, k];
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

  for (n in 1:N_obs) {
    int i = iid[n];
    int K = length_poly[i];
    vector[K] lp;
    for (k in 1:K)
      lp[k] = (k - 1) * eta[n] + d[i, k];
    log_lik[n] = lp[y[n] + 1] - log_sum_exp(lp);
  }
}
