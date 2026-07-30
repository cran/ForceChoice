// FCGDINA: forced-choice cognitive diagnostic model.
// Base CDM: DINA, DINO, ACDM, or GDINA via per-item design matrices.
// Forced-choice wrapper matches FCMIRT/FCGGUM: RANK, MOLE, and PICK are
// evaluated through the same sequential Luce/Plackett transformation.
// Structural parameters pi are jointly estimated with a Dirichlet prior.

functions {
  real fcgdina_item_eta(
    int cls, int item,
    array[] real flat_design,
    array[] int design_rows,
    array[] int design_cols,
    array[] int design_offset,
    array[] int delta_offset,
    array[,] int class_map,
    vector delta
  ) {
    int off_x = design_offset[item];
    int off_d = delta_offset[item];
    int nrow = design_rows[item];
    int ncol = design_cols[item];
    int row = class_map[item, cls];
    real eta = 0.0;

    for (c in 1:ncol) {
      eta += flat_design[off_x + (row - 1) * ncol + c] * delta[off_d + c];
    }
    return eta;
  }

  // Compute block-level pattern log-probabilities for all C classes.
  // Returns vector[C * n_obs_b] where element (c-1) * n_obs_b + t is
  // log P(observed pattern t | class c, block b).
  vector fcgdina_block_log_probs(
    // Block-specific scalars (primitive ints to allow dynamic return size)
    int K, int n_tot, int t_start, int n_obs_b,
    int obs_pat_start, int blk_item_start, int C,
    // Flat arrays for design matrices
    array[] real flat_design,
    array[] int design_rows,
    array[] int design_cols,
    array[] int design_offset,
    array[] int delta_offset,
    array[,] int class_map,
    vector delta,
    // Forced-choice data arrays
    array[] int block_items,
    array[] int patterns_total,
    array[] int obs_full_start,
    array[] int obs_full_index
  ) {
    vector[C * n_obs_b] result;

    // The normalized agree/disagree construction is algebraically a
    // Plackett-Luce model on item logits. Work on logits directly to avoid
    // repeated products and their derivatives.
    matrix[C, K] eta_mat;

    for (cls in 1:C) {
      for (k in 1:K) {
        int item = block_items[blk_item_start + k - 1];
        real eta = fcgdina_item_eta(
          cls, item, flat_design, design_rows, design_cols,
          design_offset, delta_offset, class_map, delta
        );
        eta_mat[cls, k] = eta;
      }
    }

    // Compute full-ranking log-probabilities via sequential Luce-Plackett.
    matrix[C, n_tot] lp_total;

    if (K > 1) {
      for (cls in 1:C) {
        vector[n_tot] lp_cls;

        for (t in 1:n_tot) {
          real lp = 0;
          for (pos in 1:(K - 1)) {
            int rem_len = K - pos + 1;
            vector[K] pick_lp;
            real chosen_lp = 0;

            for (h_pos in 1:rem_len) {
              int h = pos + h_pos - 1;
              int candidate = patterns_total[t_start + (t - 1) * K + h - 1];
              real term = eta_mat[cls, candidate];
              pick_lp[h_pos] = term;
              if (h_pos == 1) {
                chosen_lp = term;
              }
            }
            lp += chosen_lp - log_sum_exp(pick_lp[1:rem_len]);
          }
          lp_cls[t] = lp;
        }

        for (t in 1:n_tot) {
          lp_total[cls, t] = lp_cls[t];
        }
      }
    } else {
      // Single-item block: only one "ranking" with probability 1.
      for (cls in 1:C) {
        for (t in 1:n_tot) {
          lp_total[cls, t] = 0;
        }
      }
    }

    // Aggregate full rankings to observed patterns through precomputed lookup.
    for (cls in 1:C) {
      for (t in 1:n_obs_b) {
        int obs_idx = obs_pat_start + t - 1;
        int first_link = obs_full_start[obs_idx];
        int last_link = obs_full_start[obs_idx + 1] - 1;
        real lp_obs;

        lp_obs = lp_total[cls, obs_full_index[first_link]];
        if (last_link > first_link) {
          for (j in (first_link + 1):last_link) {
            lp_obs = log_sum_exp(lp_obs, lp_total[cls, obs_full_index[j]]);
          }
        }
        result[(cls - 1) * n_obs_b + t] = lp_obs;
      }
    }

    return result;
  }
}

data {
  int<lower=1> N;
  int<lower=1> B;
  int<lower=1> I;
  int<lower=1> D;
  int<lower=1> C;

  int<lower=1> G;
  array[G, B] int<lower=1> y_unique;
  array[G] int<lower=1> y_count;
  array[C, D] int<lower=0, upper=1> alpha_patterns;

  int<lower=1> total_design_entries;
  array[total_design_entries] real flat_design;
  array[I] int<lower=1> design_rows;
  array[I] int<lower=1> design_cols;
  array[I] int<lower=0> design_offset;

  int<lower=1> total_delta_len;
  int<lower=1> total_delta_free;
  matrix[total_delta_len, total_delta_free] delta_basis;
  array[I] int<lower=0> delta_offset;
  array[I, C] int<lower=1> class_map;

  array[B] int<lower=2> n_items;
  int<lower=1> total_block_items;
  array[total_block_items] int<lower=1, upper=I> block_items;
  array[B + 1] int<lower=1> item_start;

  array[B] int<lower=1> n_total;
  int<lower=1> total_cells_total;
  array[total_cells_total] int<lower=1> patterns_total;
  array[B + 1] int<lower=1> total_start;

  array[B] int<lower=1> n_obs;
  int<lower=1> total_obs_patterns;
  array[B + 1] int<lower=1> obs_pattern_start;
  int<lower=1> total_obs_full_links;
  array[total_obs_patterns + 1] int<lower=1> obs_full_start;
  array[total_obs_full_links] int<lower=1> obs_full_index;

  vector[total_delta_free] delta_prior_mu;
  real<lower=0> delta_prior_sigma;
  real<lower=0> pi_prior_alpha;  // Dirichlet concentration (1 = uniform)
}

transformed data {
  array[B] int<lower=1> blk_off;

  {
    int cum = 1;
    for (b in 1:B) {
      blk_off[b] = cum;
      cum += C * n_obs[b];
    }
  }
}

parameters {
  vector[total_delta_free] delta_free;
  simplex[C] pi;
}

transformed parameters {
  vector[total_delta_len] delta = delta_basis * delta_free;
  // Pre-computed pattern log-probabilities: P(observed pattern | class, block).
  // Flat storage: block-major, within each block class-major.
  // For block b, pattern t, class c:
  //   index = block_offset[b] + (c - 1) * n_obs[b] + t - 1   (0-based)
  //   or:   block_offset_1b[b] + (c - 1) * n_obs[b] + t - 1  (1-based vector)
  vector[C * total_obs_patterns] log_prob_flat;

  {
    int pos = 1;
    for (b in 1:B) {
      int n_obs_b = n_obs[b];
      vector[C * n_obs_b] block_probs = fcgdina_block_log_probs(
        // Block-specific scalars
        n_items[b], n_total[b], total_start[b], n_obs_b,
        obs_pattern_start[b], item_start[b], C,
        // Flat arrays
        flat_design, design_rows, design_cols,
        design_offset, delta_offset, class_map, delta,
        // FC data arrays
        block_items, patterns_total, obs_full_start, obs_full_index
      );
      log_prob_flat[pos : pos + C * n_obs_b - 1] = block_probs;
      pos += C * n_obs_b;
    }
  }
}

model {
  delta_free ~ normal(delta_prior_mu, delta_prior_sigma);
  pi ~ dirichlet(rep_vector(pi_prior_alpha, C));

  {
    vector[C] log_pi = log(pi);

    for (g in 1:G) {
      vector[C] class_lp = log_pi;
      for (b in 1:B) {
        int n_obs_b = n_obs[b];
        int base = blk_off[b];
        int y_gb = y_unique[g, b];
        for (cls in 1:C) {
          class_lp[cls] += log_prob_flat[base + (cls - 1) * n_obs_b + y_gb - 1];
        }
      }
      target += y_count[g] * log_sum_exp(class_lp);
    }
  }
}

generated quantities {
  matrix[C, G] class_prob_group;
  vector[G] log_lik_group;

  {
    vector[C] log_pi = log(pi);

    for (g in 1:G) {
      vector[C] lp_g = log_pi;
      for (b in 1:B) {
        int n_obs_b = n_obs[b];
        int base = blk_off[b];
        int y_gb = y_unique[g, b];
        for (cls in 1:C) {
          lp_g[cls] += log_prob_flat[base + (cls - 1) * n_obs_b + y_gb - 1];
        }
      }

      {
        vector[C] prob_g = softmax(lp_g);
        log_lik_group[g] = log_sum_exp(lp_g);
        for (cls in 1:C) {
          class_prob_group[cls, g] = prob_g[cls];
        }
      }
    }
  }
}
