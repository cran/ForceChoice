// FCMIRT model: forced-choice wrapper around a binary MIRT item model.
// model_type: 1 = 1PL, 2 = 2PL, 3 = 3PL, 4 = 4PL

functions {
  real fcmirt_obs_ll(
    int p, int blk, int y_obs,
    int model_type,
    matrix theta, matrix a, vector b_par, vector c_par, vector d_par,
    array[] int n_items, array[] int fc_type,
    array[] int block_items, array[] int item_start,
    array[] int n_total, array[] int patterns_total, array[] int total_start,
    array[] int patterns_obs, array[] int obs_start, array[] int fc_len
  ) {
    int K = n_items[blk];
    vector[K] log_agree;
    vector[K] log_disagree;

    for (k in 1:K) {
      int gi = block_items[item_start[blk] + k - 1];
      real eta = dot_product(theta[p], a[gi]) - b_par[gi];
      real p_agree;

      if (model_type == 3) {
        p_agree = c_par[gi] + (1.0 - c_par[gi]) * inv_logit(eta);
      } else if (model_type == 4) {
        p_agree = c_par[gi] + (d_par[gi] - c_par[gi]) * inv_logit(eta);
      } else {
        p_agree = inv_logit(eta);
      }

      p_agree = fmin(1.0 - 1e-12, fmax(1e-12, p_agree));
      log_agree[k] = log(p_agree);
      log_disagree[k] = log1m(p_agree);
    }

    {
      int n_tot = n_total[blk];
      int t_start = total_start[blk];
      vector[n_tot] lp_total;

      for (t in 1:n_tot) {
        real lp = 0;
        if (K > 1) {
          for (pos in 1:(K - 1)) {
            int rem_len = K - pos + 1;
            vector[K] pick_lp;
            real chosen_lp = 0;

            for (h_pos in 1:rem_len) {
              int h = pos + h_pos - 1;
              int candidate = patterns_total[t_start + (t - 1) * K + h - 1];
              real term = log_agree[candidate];

              for (g in pos:K) {
                if (g != h) {
                  int other = patterns_total[t_start + (t - 1) * K + g - 1];
                  term += log_disagree[other];
                }
              }

              pick_lp[h_pos] = term;
              if (h_pos == 1) {
                chosen_lp = term;
              }
            }

            lp += chosen_lp - log_sum_exp(pick_lp[1:rem_len]);
          }
        }
        lp_total[t] = lp;
      }

      {
        int fc_t = fc_type[blk];
        int o_start = obs_start[blk];
        int fcl = fc_len[blk];
        real lp_obs = negative_infinity();

        if (fc_t == 1) {
          return fmax(lp_total[y_obs], log(1e-16));
        }

        for (t in 1:n_tot) {
          int consistent = 0;
          if (fc_t == 2) {
            int first_t = patterns_total[t_start + (t - 1) * K];
            int last_t  = patterns_total[t_start + (t - 1) * K + K - 1];
            int first_o = patterns_obs[o_start + (y_obs - 1) * fcl];
            int last_o  = patterns_obs[o_start + (y_obs - 1) * fcl + 1];
            consistent = (first_t == first_o) && (last_t == last_o);
          } else {
            int first_t = patterns_total[t_start + (t - 1) * K];
            int first_o = patterns_obs[o_start + (y_obs - 1) * fcl];
            consistent = (first_t == first_o);
          }
          if (consistent) {
            lp_obs = log_sum_exp(lp_obs, lp_total[t]);
          }
        }

        return fmax(lp_obs, log(1e-16));
      }
    }
  }
}

data {
  int<lower=1> N_obs;
  int<lower=1> N;
  int<lower=1> B;
  int<lower=1> I;
  int<lower=1> D;
  int<lower=1, upper=4> model_type;

  array[N_obs] int<lower=1, upper=N> pid;
  array[N_obs] int<lower=1, upper=B> bid;
  array[N_obs] int<lower=1> y;

  matrix[I, D] Q_matrix;

  array[B] int<lower=1> n_items;
  array[B] int<lower=1, upper=3> fc_type;
  int<lower=1> total_block_items;
  array[total_block_items] int<lower=1, upper=I> block_items;
  array[B + 1] int<lower=1> item_start;

  array[B] int<lower=1> n_total;
  int<lower=1> total_cells_total;
  array[total_cells_total] int<lower=1> patterns_total;
  array[B + 1] int<lower=1> total_start;

  array[B] int<lower=1> n_obs;
  int<lower=1> total_cells_obs;
  array[total_cells_obs] int<lower=1> patterns_obs;
  array[B + 1] int<lower=1> obs_start;
  array[B] int<lower=1> fc_len;

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
  int b_free_size = I - B;
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
  if (model_type >= 3) {
    c_size = I;
  }
  if (model_type == 4) {
    d_size = I;
  }
}

parameters {
  vector<lower=0>[a_free_size] a_free;
  vector[b_free_size] b_free;
  vector<lower=c_mu, upper=c_sigma>[c_size] c;
  vector<lower=d_mu, upper=d_sigma>[d_size] d;
  cholesky_factor_corr[D > 1 ? D : 2] L_Corr;
  array[N] vector[D] theta_raw;
}

transformed parameters {
  matrix<lower=0>[I, D] a;
  vector[I] b;
  vector[I] c_full;
  vector[I] d_full;
  matrix[N, D] theta;

  if (model_type == 1) {
    for (i in 1:I) {
      for (j in 1:D) {
        a[i, j] = Q_matrix[i, j];
      }
    }
  } else {
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

  {
    int idx = 1;
    for (blk in 1:B) {
      int K = n_items[blk];
      real b_sum = 0.0;
      for (k in 1:(K - 1)) {
        int gi = block_items[item_start[blk] + k - 1];
        b[gi] = b_free[idx];
        b_sum += b_free[idx];
        idx += 1;
      }
      b[block_items[item_start[blk] + K - 1]] = -b_sum;
    }
  }

  if (model_type >= 3) {
    c_full = c;
  } else {
    c_full = rep_vector(0.0, I);
  }

  if (model_type == 4) {
    d_full = d;
  } else {
    d_full = rep_vector(1.0, I);
  }

  for (n in 1:N) {
    if (D > 1) {
      theta[n] = (theta_mu + L_Corr * theta_raw[n])';
    } else {
      theta[n, 1] = theta_mu[1] + theta_raw[n][1];
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

  for (obs in 1:N_obs) {
    target += fcmirt_obs_ll(
      pid[obs], bid[obs], y[obs], model_type,
      theta, a, b, c_full, d_full,
      n_items, fc_type,
      block_items, item_start,
      n_total, patterns_total, total_start,
      patterns_obs, obs_start, fc_len
    );
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

  for (obs in 1:N_obs) {
    log_lik[obs] = fcmirt_obs_ll(
      pid[obs], bid[obs], y[obs], model_type,
      theta, a, b, c_full, d_full,
      n_items, fc_type,
      block_items, item_start,
      n_total, patterns_total, total_start,
      patterns_obs, obs_start, fc_len
    );
  }
}
