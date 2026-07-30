// FCDCM: forced-choice diagnostic classification model.
// Attribute profiles are marginalized exactly conditional on a higher-order
// latent theta. Equal-condensation probabilities are fixed at 0.5.

functions {
  matrix fcdcm_response_logprob(int C, int N, int B, array[] int block_case,
                                array[] int y, vector eta0, vector etaAB) {
    matrix[C, N] logprob = rep_matrix(0.0, C, N);
    vector[B] log_eta0;
    vector[B] log_eta_high;
    vector[B] log_not_eta0;
    vector[B] log_not_eta_high;
    real log_half = log(0.5);

    for (b in 1:B) {
      real eta0_safe = fmin(0.5 - 1e-12, fmax(1e-12, eta0[b]));
      real eta_high_safe = fmin(1.0 - 1e-12, fmax(1e-12, 0.5 + etaAB[b]));
      log_eta0[b] = log(eta0_safe);
      log_eta_high[b] = log(eta_high_safe);
      log_not_eta0[b] = log1m(eta0_safe);
      log_not_eta_high[b] = log1m(eta_high_safe);
    }

    for (cls in 1:C) {
      int class_offset = (cls - 1) * B;
      for (b in 1:B) {
        int case_id = block_case[class_offset + b];
        int y_offset = (b - 1) * N;
        real logprob_y0;
        real logprob_y1;

        if (case_id == 0) {
          // zeta_a < zeta_b: choose A with eta0
          logprob_y0 = log_not_eta0[b];
          logprob_y1 = log_eta0[b];
        } else if (case_id == 1) {
          // zeta_a == zeta_b == 0: fixed at 0.5
          logprob_y0 = log_half;
          logprob_y1 = log_half;
        } else if (case_id == 2) {
          // zeta_a > zeta_b: choose A with 0.5 + etaAB
          logprob_y0 = log_not_eta_high[b];
          logprob_y1 = log_eta_high[b];
        } else {
          // zeta_a == zeta_b == 1: fixed at 0.5
          logprob_y0 = log_half;
          logprob_y1 = log_half;
        }

        for (n in 1:N) {
          if (y[y_offset + n] == 1) {
            logprob[cls, n] += logprob_y1;
          } else {
            logprob[cls, n] += logprob_y0;
          }
        }
      }
    }

    return logprob;
  }
}

data {
  int<lower=1> N;
  int<lower=1> B;
  int<lower=1> I;
  int<lower=1> D;
  int<lower=1> C;

  array[N * B] int<lower=0, upper=1> y;

  array[C * D] int<lower=0, upper=1> alpha_patterns;
  array[C * I] int<lower=0, upper=1> zeta_patterns;
  array[B * 2] int<lower=1, upper=I> patterns;

  real delta1_mu;
  real<lower=0> delta1_sigma;
  real delta0_mu;
  real<lower=0> delta0_sigma;
  real<lower=0> eta0_mu;
  real<lower=0> eta0_sigma;
  real<lower=0> etaAB_mu;
  real<lower=0> etaAB_sigma;
}

transformed data {
  array[C * B] int<lower=0, upper=3> block_case;
  matrix[C, D] alpha_matrix;
  matrix[C, D] alpha_nonmatrix;

  for (cls in 1:C) {
    int class_offset = (cls - 1) * I;
    int block_offset = (cls - 1) * B;

    for (d in 1:D) {
      alpha_matrix[cls, d] = alpha_patterns[(cls - 1) * D + d];
      alpha_nonmatrix[cls, d] = 1 - alpha_matrix[cls, d];
    }

    for (b in 1:B) {
      int item_a = patterns[(b - 1) * 2 + 1];
      int item_b = patterns[(b - 1) * 2 + 2];
      int zeta_a = zeta_patterns[class_offset + item_a];
      int zeta_b = zeta_patterns[class_offset + item_b];

      if (zeta_a < zeta_b) {
        block_case[block_offset + b] = 0;
      } else if (zeta_a > zeta_b) {
        block_case[block_offset + b] = 2;
      } else if (zeta_a == 0) {
        block_case[block_offset + b] = 1;
      } else {
        block_case[block_offset + b] = 3;
      }
    }
  }
}

parameters {
  vector[N] theta;
  vector<lower=0>[D] delta1;
  vector[D] delta0;
  vector<lower=0, upper=0.5>[B] eta0;
  vector<lower=0, upper=0.5>[B] etaAB;
}

model {
  matrix[C, N] response_logprob = fcdcm_response_logprob(C, N, B,
                                                         block_case, y,
                                                         eta0, etaAB);

  theta ~ std_normal();
  delta1 ~ lognormal(delta1_mu, delta1_sigma);
  delta0 ~ normal(delta0_mu, delta0_sigma);
  eta0 ~ normal(eta0_mu, eta0_sigma);
  etaAB ~ normal(etaAB_mu, etaAB_sigma);

  for (n in 1:N) {
    vector[D] eta = delta1 .* (rep_vector(theta[n], D) - delta0);
    vector[D] log_mastery = log_inv_logit(eta);
    vector[D] log_nonmastery = log_inv_logit(-eta);
    vector[C] alpha_lp = alpha_matrix * log_mastery +
      alpha_nonmatrix * log_nonmastery;
    target += log_sum_exp(alpha_lp + col(response_logprob, n));
  }
}

generated quantities {
  vector[N] log_lik;

  {
    matrix[C, N] response_logprob = fcdcm_response_logprob(C, N, B,
                                                           block_case, y,
                                                           eta0, etaAB);

    for (n in 1:N) {
      vector[D] eta = delta1 .* (rep_vector(theta[n], D) - delta0);
      vector[D] log_mastery = log_inv_logit(eta);
      vector[D] log_nonmastery = log_inv_logit(-eta);
      vector[C] alpha_lp = alpha_matrix * log_mastery +
        alpha_nonmatrix * log_nonmastery;
      log_lik[n] = log_sum_exp(alpha_lp + col(response_logprob, n));
    }
  }
}
