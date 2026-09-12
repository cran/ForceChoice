// FCGGUM Model — Forced-Choice Generalized Graded Unfolding Model
//
// This is a forced-choice wrapper around the MGGUM.  At the item level the
// probability of "agreeing with" (endorsing) item i is obtained from the
// MGGUM category-2 probability (MUPP-style, C=K-1, M=2C+1=2K-1).
// Full rankings are modeled as sequential Luce/PICK decisions over the
// remaining statements in each block.
//
//   RANK  – full ranking of all items within a block
//   MOLE  – most-least (best-worst) partial ranking
//   PICK  – single-best selection
//
// Identification:
//   - theta_raw ~ std_normal() with theta_mu = 0 (location + scale of theta)
//   - delta side constraints (Q=1: delta >= 0, Q=-1: delta <= 0);
//     Q sign does not make discrimination negative.
//   - tau ordered: all items are binary (Ki=2), so a single positive neg_tau
//     per item (tau[i,2] = -neg_tau[i] < 0, no intra-item ordering needed)

functions {
  // Return the observed-choice log probability (pure function, no target()).
  real fcggum_obs_ll(
    int p, int b, int y_obs,
    matrix theta, matrix a, matrix delta, matrix tau,
    int D,
    array[] int n_items, array[] int fc_type,
    array[] int block_items, array[] int item_start,
    array[] int length_poly,
    array[] int n_total, array[] int patterns_total, array[] int total_start,
    array[] int patterns_obs, array[] int obs_start, array[] int fc_len
  ) {
    int K = n_items[b];
    vector[K] log_agree;
    vector[K] log_disagree;

    for (k in 1:K) {
      int gi = block_items[item_start[b] + k - 1];
      int Ki = length_poly[gi];
      int Mi = 2 * Ki - 1;
      real sum_a_i = sum(a[gi]);

      real r_sq = 0;
      for (d in 1:D)
        r_sq += square(a[gi, d]) * square(theta[p, d] - delta[gi, d]);
      real r = sqrt(r_sq + 1e-12);

      real cum_psi = 0;
      vector[Ki] lp_cat;
      for (cat in 1:Ki) {
        cum_psi += tau[gi, cat] * sum_a_i;
        lp_cat[cat] = log_sum_exp(
          (cat - 1) * r - cum_psi,
          (Mi - cat + 1) * r - cum_psi
        );
      }
      {
        real p_raw = exp(lp_cat[2] - log_sum_exp(lp_cat));
        real p_agree = fmin(1.0 - 1e-12, fmax(1e-12, p_raw));
        log_agree[k] = log(p_agree);
        log_disagree[k] = log1m(p_agree);
      }
    }

    int n_tot = n_total[b];
    vector[n_tot] lp_total;
    int t_start = total_start[b];

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
            if (h_pos == 1)
              chosen_lp = term;
          }

          lp += chosen_lp - log_sum_exp(pick_lp[1:rem_len]);
        }
      }
      lp_total[t] = lp;
    }

    int fc_t = fc_type[b];
    int o_start = obs_start[b];
    int fcl = fc_len[b];
    real lp_obs = negative_infinity();

    if (fc_t == 1) {
      return fmax(lp_total[y_obs], log(1e-16));
    } else {
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
        if (consistent)
          lp_obs = log_sum_exp(lp_obs, lp_total[t]);
      }

      return fmax(lp_obs, log(1e-16));
    }
  }

  // Convenience wrapper that increments target() with the result of fcggum_obs_ll.
  void fcggum_obs_lp(
    int p, int b, int y_obs,
    matrix theta, matrix a, matrix delta, matrix tau,
    int D,
    array[] int n_items, array[] int fc_type,
    array[] int block_items, array[] int item_start,
    array[] int length_poly,
    array[] int n_total, array[] int patterns_total, array[] int total_start,
    array[] int patterns_obs, array[] int obs_start, array[] int fc_len
  ) {
    target += fcggum_obs_ll(
      p, b, y_obs, theta, a, delta, tau,
      D, n_items, fc_type, block_items, item_start,
      length_poly, n_total, patterns_total, total_start,
      patterns_obs, obs_start, fc_len
    );
  }
}

data {
  // ---- Observation dimensions ----
  int<lower=1> N_obs;              // N * B
  int<lower=1> N;
  int<lower=1> B;
  int<lower=1> I;
  int<lower=1> D;
  int<lower=1> max_poly;           // = 2 for FCGGUM

  // ---- Observations (long format) ----
  array[N_obs] int<lower=1, upper=N> pid;
  array[N_obs] int<lower=1, upper=B> bid;
  array[N_obs] int<lower=1>         y;

  // ---- Item metadata ----
  array[I] int<lower=2, upper=max_poly> length_poly;   // = 2 per item
  matrix[I, D] Q_matrix;                               // -1/1 = active delta side; 0 = inactive

  // ---- Block structure ----
  array[B] int<lower=1>     n_items;
  array[B] int<lower=1, upper=3> fc_type;    // 1 = RANK, 2 = MOLE, 3 = PICK
  int<lower=1>               total_block_items;
  array[total_block_items] int<lower=1, upper=I> block_items;
  array[B + 1] int<lower=1>  item_start;

  // ---- Full ranking patterns (for sequential ranking probabilities) ----
  array[B] int<lower=1>     n_total;
  int<lower=1>               total_cells_total;
  array[total_cells_total] int<lower=1> patterns_total;
  array[B + 1] int<lower=1>  total_start;

  // ---- FC-type patterns (for observation) ----
  array[B] int<lower=1>     n_obs;
  int<lower=1>               total_cells_obs;
  array[total_cells_obs] int<lower=1> patterns_obs;
  array[B + 1] int<lower=1>  obs_start;
  array[B] int<lower=1>     fc_len;

  // ---- Prior hyperparameters (overridable via control) ----
  real a_mu;                   // lognormal location for a (default 0 -> median=1)
  real<lower=0> a_sigma;       // lognormal scale    for a (default 0.5)
  real delta_pos_mu;           // normal location for positive delta (default 1)
  real<lower=0> delta_pos_sigma; // normal scale   for positive delta (default 0.5)
  real delta_neg_mu;           // normal location for negative delta (default -1)
  real<lower=0> delta_neg_sigma; // normal scale   for negative delta (default 0.5)
  real tau_mu;                 // lognormal location for neg_tau (default 0 -> median=1)
  real<lower=0> tau_sigma;     // lognormal scale    for neg_tau (default 0.5)

  vector[D] theta_mu;
}

transformed data {
  int a_free_size = 0;
  int delta_pos_size = 0;
  int delta_neg_size = 0;
  for (i in 1:I) {
    for (j in 1:D) {
      if (Q_matrix[i, j] != 0) {
        a_free_size += 1;
        if (Q_matrix[i, j] == 1)
          delta_pos_size += 1;
        else
          delta_neg_size += 1;
      }
    }
  }
}

parameters {
  vector<lower=0>[a_free_size] a_free;          // positivity constraint
  vector<upper=0>[delta_neg_size] delta_neg;     // Q=-1: negative-side delta
  vector<lower=0>[delta_pos_size] delta_pos;     // Q=1: positive-side delta
  vector<lower=0>[I] neg_tau;                    // one positive scalar per item (Ki=2)
  cholesky_factor_corr[D] L_Corr;
  array[N] vector[D] theta_raw;
}

transformed parameters {
  matrix<lower=0>[I, D] a;     // lower=0 allows zero entries for Q=0
  matrix[I, D] delta;
  matrix[I, max_poly] tau;
  matrix[N, D] theta;

  // ---- Fill a matrix ----
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

  // ---- Fill delta matrix ----
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

  // ---- Fill tau matrix ----
  // All items are binary (Ki=2): tau[i,1] = 0, tau[i,2] = -neg_tau[i] < 0
  for (i in 1:I) {
    tau[i, 1] = 0.0;
    tau[i, 2] = -neg_tau[i];
    if (max_poly > 2) {
      for (k in 3:max_poly)
        tau[i, k] = 0.0;
    }
  }

  // ---- Compute theta ----
  for (n in 1:N) {
    if (D > 1)
      theta[n] = (theta_mu + L_Corr * theta_raw[n])';
    else
      theta[n, 1] = theta_mu[1] + theta_raw[n][1];
  }
}

model {
  // ── Priors (smooth, HMC-friendly; hyperparameters overridable via control) ──
  a_free    ~ lognormal(a_mu, a_sigma);
  delta_neg ~ normal(delta_neg_mu, delta_neg_sigma);
  delta_pos ~ normal(delta_pos_mu, delta_pos_sigma);
  neg_tau   ~ lognormal(tau_mu, tau_sigma);

  if (D > 1)
    L_Corr ~ lkj_corr_cholesky(1);

  for (n in 1:N)
    theta_raw[n] ~ std_normal();

  // ── Likelihood ──
  // Forced-choice GGUM: item-level endorsement probs from MGGUM (category 2 = "agree"),
  // combined into ranking-pattern probabilities via sequential Luce/PICK decisions.
  for (obs in 1:N_obs) {
    fcggum_obs_lp(
      pid[obs], bid[obs], y[obs],
      theta, a, delta, tau,
      D,
      n_items, fc_type,
      block_items, item_start,
      length_poly,
      n_total, patterns_total, total_start,
      patterns_obs, obs_start, fc_len
    );
  }
}

generated quantities {
  matrix[D, D] Corr;
  vector[N_obs] log_lik;

  if (D > 1)
    Corr = multiply_lower_tri_self_transpose(L_Corr);
  else
    Corr[1, 1] = 1.0;

  for (obs in 1:N_obs) {
    log_lik[obs] = fcggum_obs_ll(
      pid[obs], bid[obs], y[obs],
      theta, a, delta, tau,
      D,
      n_items, fc_type,
      block_items, item_start,
      length_poly,
      n_total, patterns_total, total_start,
      patterns_obs, obs_start, fc_len
    );
  }
}
