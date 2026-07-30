// TIRT Model — Thurstonian IRT for Forced-Choice Data
//
// Pairwise probit comparison of latent utility differences. For each
// forced-choice pair, a univariate normal latent difference is formed from
// the statement-level factor loadings (lambda) and person thetas; the binary
// pairwise outcome follows a probit link.
//
// Supports RANK, MOLE, and PICK forced-choice formats via sequential
// Luce/Plackett transformation of the pairwise probabilities.
//
// Identification:
//   - theta_raw ~ std_normal() with theta_mu = 0 (location + scale of theta)
//   - Corr ~ lkj_corr_cholesky(1)
//   - Psi2 (uniquenesses) bounded in (0, 1)
data {
  int<lower=1> N_obs;
  int<lower=1> N_pair;
  int<lower=1> N_person;
  int<lower=2> I_states;
  int<lower=2> D;

  array[N_obs] int<lower=0, upper=1> response;
  array[N_obs] int<lower=1, upper=N_person> Idx_person;
  array[N_obs] int<lower=1, upper=I_states> Idx_item_i;
  array[N_obs] int<lower=1, upper=I_states> Idx_item_k;
  array[N_obs] int<lower=1, upper=N_pair> Idx_pair;
  matrix[I_states, D] Q_matrix;

  real psi_mu;
  real<lower=0> psi_sigma;
  real gamma_mu;
  real<lower=0> gamma_sigma;
  vector[D] theta_mu;

  real<lower=0> lambda_alpha;
  real<lower=0> lambda_beta;

  int<lower=0> N_lambda_est;
  array[N_lambda_est] int<lower=1, upper=I_states> Idx_lambda_est;
  int<lower=0> N_lambda_fix;
  array[N_lambda_fix] int<lower=1, upper=I_states> Idx_lambda_fix;
  vector[N_lambda_fix] lambda_fix_val;

  int<lower=0> N_psi_fix;
  array[N_psi_fix] int<lower=1, upper=I_states> Idx_psi_fix;
  vector<lower=0>[N_psi_fix] psi_fix_val;
  int<lower=0> N_psi_est;
  array[N_psi_est] int<lower=1, upper=I_states> Idx_psi_est;
  int<lower=0> N_psi_equal;
  array[N_psi_equal] int<lower=1, upper=I_states> Idx_psi_equal;
  array[N_psi_equal] int<lower=1, upper=I_states> Idx_psi_orig;

  int<lower=0> N_gamma_fix;
  array[N_gamma_fix] int<lower=1, upper=N_pair> Idx_gamma_fix;
  vector[N_gamma_fix] gamma_fix_val;
  int<lower=0> N_gamma_est;
  array[N_gamma_est] int<lower=1, upper=N_pair> Idx_gamma_est;
}

parameters {
  matrix[D, N_person] z_theta;
  cholesky_factor_corr[D] L_theta;
  vector<lower=lambda_alpha, upper=lambda_beta>[N_lambda_est] lambda_est;
  vector<lower=0>[N_psi_est] psi_est;
  vector[N_gamma_est] gamma_est;
}

transformed parameters {
  matrix[N_person, D] theta;
  vector<lower=0>[I_states] lambda;
  vector<lower=0>[I_states] psi;

  theta = (L_theta * z_theta)';

  lambda[Idx_lambda_est] = lambda_est;
  lambda[Idx_lambda_fix] = lambda_fix_val;

  psi[Idx_psi_fix] = psi_fix_val;
  psi[Idx_psi_est] = psi_est;
  psi[Idx_psi_equal] = psi[Idx_psi_orig];

  vector[N_pair] gamma;
  gamma[Idx_gamma_fix] = gamma_fix_val;
  gamma[Idx_gamma_est] = gamma_est;
}

model {
  to_vector(z_theta) ~ std_normal();
  L_theta ~ lkj_corr_cholesky(1);
  lambda_est ~ uniform(lambda_alpha, lambda_beta);
  psi_est ~ normal(psi_mu, psi_sigma);
  gamma_est ~ normal(gamma_mu, gamma_sigma);

  vector[N_obs] mu_diff;
  vector[N_obs] sum_psi;

  for (n in 1:N_obs) {
    int i = Idx_item_i[n];
    int k = Idx_item_k[n];
    int p = Idx_person[n];

    mu_diff[n] = lambda[i] * dot_product(Q_matrix[i], theta[p]) -
                 lambda[k] * dot_product(Q_matrix[k], theta[p]);

    sum_psi[n] = sqrt(psi[i]^2 + psi[k]^2 + 1e-12);
  }

  for (n in 1:N_obs) {
    real z = (mu_diff[n] - gamma[Idx_pair[n]]) / sum_psi[n];
    if (response[n] == 1) {
      target += normal_lcdf(z | 0, 1);
    } else {
      target += normal_lccdf(z | 0, 1);
    }
  }
}

generated quantities {
  corr_matrix[D] Corr;
  vector[N_obs] log_lik;

  Corr = multiply_lower_tri_self_transpose(L_theta);

  for (n in 1:N_obs) {
    real mu_diff = lambda[Idx_item_i[n]] * dot_product(Q_matrix[Idx_item_i[n]], theta[Idx_person[n]]) -
                   lambda[Idx_item_k[n]] * dot_product(Q_matrix[Idx_item_k[n]], theta[Idx_person[n]]);
    real sum_psi = sqrt(psi[Idx_item_i[n]]^2 + psi[Idx_item_k[n]]^2 + 1e-12);
    real z = (mu_diff - gamma[Idx_pair[n]]) / sum_psi;
    if (response[n] == 1) {
      log_lik[n] = normal_lcdf(z | 0, 1);
    } else {
      log_lik[n] = normal_lccdf(z | 0, 1);
    }
  }
}
