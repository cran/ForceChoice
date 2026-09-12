// Unified MIRT model (m1pl / m2pl / m3pl / m4pl)
// model_type: 1 = 1PL, 2 = 2PL, 3 = 3PL, 4 = 4PL
data {
  int<lower=1> N_obs;
  int<lower=1> N;
  int<lower=1> I;
  int<lower=1> D;
  int<lower=1, upper=4> model_type;

  array[N_obs] int<lower=0, upper=1> y;
  array[N_obs] int<lower=1, upper=N> pid;
  array[N_obs] int<lower=1, upper=I> iid;

  matrix[I, D] Q_matrix;

  real a_mu;
  real<lower=0> a_sigma;
  real b_mu;
  real<lower=0> b_sigma;
  real c_mu;
  real<lower=0> c_sigma;
  real d_mu;
  real<lower=0> d_sigma;

  vector[D] theta_mu;
}

transformed data {
  int a_free_size = 0;
  int c_size = 0;
  int d_size = 0;

  if (model_type >= 2) {
    for (i in 1:I) {
      for (j in 1:D) {
        if (Q_matrix[i, j] == 1) {
          a_free_size += 1;
        }
      }
    }
  }
  if (model_type >= 3) c_size = I;
  if (model_type == 4) d_size = I;
}

parameters {
  vector<lower=0>[a_free_size] a_free;
  vector<lower=-4, upper=4>[I] b;
  vector<lower=c_mu, upper=c_sigma>[c_size] c;
  vector<lower=d_mu, upper=d_sigma>[d_size] d;
  cholesky_factor_corr[D > 1 ? D : 2] L_Corr;
  array[N] vector[D] theta_raw;
}

transformed parameters {
  matrix[I, D] a;
  matrix[N, D] theta;
  vector[N_obs] eta;
  vector[N_obs] p;

  // --- build a matrix ---
  if (model_type == 1) {
    for (i in 1:I) {
      for (j in 1:D) {
        a[i, j] = 1.0;
      }
    }
  } else {
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
  }

  // --- build theta ---
  for (n in 1:N) {
    if (D > 1) {
      theta[n] = (theta_mu + L_Corr * theta_raw[n])';
    } else {
      theta[n, 1] = theta_mu[1] + theta_raw[n][1];
    }
  }

  // --- build eta ---
  for (n in 1:N_obs) {
    eta[n] = dot_product(theta[pid[n]], a[iid[n]]) - b[iid[n]];
  }

  // --- build p ---
  if (model_type == 3) {
    for (n in 1:N_obs) {
      real p_raw = c[iid[n]] + (1.0 - c[iid[n]]) * inv_logit(eta[n]);
      p[n] = fmin(1.0 - 1e-12, fmax(1e-12, p_raw));
    }
  } else if (model_type == 4) {
    for (n in 1:N_obs) {
      real p_raw = c[iid[n]] + (d[iid[n]] - c[iid[n]]) * inv_logit(eta[n]);
      p[n] = fmin(1.0 - 1e-12, fmax(1e-12, p_raw));
    }
  } else {
    for (n in 1:N_obs) {
      p[n] = fmin(1.0 - 1e-12, fmax(1e-12, inv_logit(eta[n])));
    }
  }
}

model {
  if (model_type >= 2) {
    a_free ~ lognormal(a_mu, a_sigma);
  }
  b ~ normal(b_mu, b_sigma);
  if (model_type >= 3) {
    c ~ uniform(c_mu, c_sigma);
  }
  if (model_type == 4) {
    d ~ uniform(d_mu, d_sigma);
  }

  if (D > 1) {
    L_Corr ~ lkj_corr_cholesky(1);
  }

  for (n in 1:N) {
    theta_raw[n] ~ std_normal();
  }

  if (model_type <= 2) {
    y ~ bernoulli_logit(eta);
  } else {
    y ~ bernoulli(p);
  }
}

generated quantities {
  matrix[D, D] Corr;
  vector[N_obs] log_lik;

  if (D > 1) {
    Corr = multiply_lower_tri_self_transpose(L_Corr);
  } else {
    Corr[1, 1] = 1.0;
  }

  for (n in 1:N_obs) {
    if (model_type <= 2) {
      log_lik[n] = bernoulli_logit_lpmf(y[n] | eta[n]);
    } else {
      log_lik[n] = bernoulli_lpmf(y[n] | p[n]);
    }
  }
}
