#' Fit TIRT Model via Stan
#'
#' @param data An N x B data frame where each column is a block and each cell
#'   contains the ranking order as a character string (e.g., "3>1>2").
#' @param Q.matrix An I x D matrix with values -1, 0, or 1. Each row has exactly
#'   one nonzero element: 1 = positive loading, -1 = reverse scoring.
#' @param block.items A list of length B, each element an integer vector of
#'   item indices belonging to that block.
#' @param fc.type Character vector of forced-choice types: \code{"RANK"} (default),
#'   \code{"MOLE"}, or \code{"PICK"}.
#' @param control.model Optional named list of model-level controls and
#'   hyperparameters. Supported entries include \code{lambda.alpha},
#'   \code{lambda.beta}, \code{psi.mu}, \code{psi.sigma}, \code{gamma.mu},
#'   \code{gamma.sigma}, and \code{theta.mu}.
#' @param control.method Optional named list of method-specific controls. For
#'   Stan this includes \code{cores}, \code{vis}, \code{seed}, \code{chains},
#'   \code{iter}, \code{warmup}, \code{thin}, \code{init}, \code{algorithm},
#'   and Stan sampler controls.
#'
#' @noRd
fit.TIRT.stan <- function(data, Q.matrix, block.items = NULL, fc.type = "RANK",
                     control.model = NULL,
                     control.method = NULL,
                       .call = NULL) {
  call <- if (is.null(.call)) match.call() else .call
  control <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  stan.method <- fc_stan_method_control(control.method)
  chains <- stan.method$chains
  iter <- stan.method$iter
  warmup <- stan.method$warmup
  thin <- stan.method$thin
  init <- stan.method$init
  algorithm <- stan.method$algorithm

  mcmc_env <- setup_mcmc_env(common.method$cores)
  on.exit(restore_mcmc_env(mcmc_env), add = TRUE)

  stan_control <- build_stan_control(algorithm, control.method)

  N.block.guess <- if (is.list(block.items)) length(block.items) else ncol(data)
  if (length(fc.type) != N.block.guess) fc.type <- rep(fc.type[1], N.block.guess)

  if (is.null(block.items)) block.items <- get.block.items.from.data(data)

  result <- get.response.from.data.TIRT(data, block.items, fc.type = fc.type)
  response <- result$response
  pairs.value <- result$pairs.value

  if (any(is.na(response))) {
    first_missing <- which(is.na(response), arr.ind = TRUE)[1L, ]
    stop(
      sprintf(
        "TIRT responses must not contain missing values (row %d, column %d).",
        first_missing[1L], first_missing[2L]
      ),
      call. = FALSE
    )
  }
  if (!all(response %in% c(0, 1))) {
    stop("TIRT responses must contain only 0 and 1.", call. = FALSE)
  }

  if (!is.matrix(Q.matrix)) {
    stop("'Q.matrix' must be a matrix.", call. = FALSE)
  }
  if (!all(Q.matrix %in% c(-1, 0, 1))) {
    stop("'Q.matrix' must contain only -1, 0, or 1.", call. = FALSE)
  }

  I.states <- nrow(Q.matrix)
  D <- ncol(Q.matrix)
  if (I.states < 1) {
    stop("'Q.matrix' must have at least one row (item).", call. = FALSE)
  }
  if (D < 2) {
    stop("At least two trait dimensions must be specified in Q.matrix (TIRT requires D >= 2).", call. = FALSE)
  }
  if (any(rowSums(abs(Q.matrix)) < 1)) {
    stop("Each item must measure at least one trait (row sums of abs('Q.matrix') must be >= 1).", call. = FALSE)
  }

  Idx_lambda_est_raw <- which(rowSums(abs(Q.matrix)) > 0)

  sm <- stanmodels[["TIRT"]]
  if (is.null(sm)) stop("Stan model 'TIRT' not found in the package namespace!", call. = FALSE)

  N.person <- nrow(response)
  if (is.null(N.person) || N.person < 1) {
    stop("Response matrix must have at least one row (person).", call. = FALSE)
  }

  if (!is.list(block.items)) {
    stop("'block.items' must be a list.", call. = FALSE)
  }
  if (length(block.items) < 1) {
    stop("'block.items' must contain at least one block.", call. = FALSE)
  }
  if (!all(sapply(block.items, function(x) length(x) >= 2))) {
    stop("Each block in 'block.items' must contain at least 2 items to form pairs.", call. = FALSE)
  }
  if (!all(sapply(block.items, is.integer) | sapply(block.items, is.numeric))) {
    stop("Elements inside 'block.items' must be integer or numeric item indices.", call. = FALSE)
  }
  all_items_in_blocks <- unlist(block.items)
  if (!all(all_items_in_blocks %in% 1:I.states)) {
    stop("Item indices in 'block.items' must be within the range of 1 to I.states (number of items in Q.matrix).", call. = FALSE)
  }
  if (length(all_items_in_blocks) != length(unique(all_items_in_blocks))) {
    stop("Item indices in 'block.items' must be unique across all blocks. An item cannot appear in multiple blocks.", call. = FALSE)
  }
  if (length(all_items_in_blocks) != I.states) {
    stop("The total number of items in 'block.items' must exactly match the number of items in Q.matrix (I.states).", call. = FALSE)
  }

  I.block  <- length(block.items[[1]])
  if (!all(sapply(block.items, length) == I.block)) {
    stop("All blocks in 'block.items' must have the same number of items.", call. = FALSE)
  }
  N.block    <- length(block.items)

  # Pre-compute per-block pair counts and total
  I.block.vec  <- sapply(block.items, length)
  block.pair.counts <- integer(N.block)
  for (b in 1:N.block) {
    k <- I.block.vec[b]
    if (fc.type[b] == "RANK") {
      block.pair.counts[b] <- choose(k, 2)
    } else if (fc.type[b] == "MOLE") {
      block.pair.counts[b] <- 2L * k - 3L
    } else if (fc.type[b] == "PICK") {
      block.pair.counts[b] <- k - 1L
    }
  }
  N.pairs.theoretical <- sum(block.pair.counts)

  N.pairs      <- ncol(response)
  if (N.pairs != N.pairs.theoretical) {
    stop(sprintf(
      "Column dimension of response matrix (%d) does not match the theoretical number of pairs derived from block.items and fc.type (%d). Please ensure response is in wide format (N.person x N.pairs) matrix.",
      N.pairs, N.pairs.theoretical
    ), call. = FALSE)
  }

  # Build pairs.matrix of ALL possible pairs (for gamma indexing)
  N.pairs.all <- sum(choose(I.block.vec, 2))
  pairs.matrix <- matrix(0L, N.pairs.all, 2)
  idx <- 0L
  for (b in 1:N.block) {
    items.b   <- block.items[[b]]
    pirs.cur  <- t(combn(items.b, 2))
    pairs.cur.length <- nrow(pirs.cur)
    pairs.matrix[(idx+1):(idx + pairs.cur.length), ] <- pirs.cur
    idx <- idx + pairs.cur.length
  }

  # Fast pair-to-gamma-index lookup (string key -> row in pairs.matrix)
  pair.keys <- paste0(pairs.matrix[,1L], "_", pairs.matrix[,2L])

  N.obs  <- N.person * N.pairs

  Y          <- integer(N.obs)
  Idx_person <- integer(N.obs)
  Idx_item_i <- integer(N.obs)
  Idx_item_k <- integer(N.obs)
  Idx_pair   <- integer(N.obs)

  idx <- 0L
  col_offset <- 0L
  for (b in 1:N.block) {
    n_pairs.b <- block.pair.counts[b]
    if (fc.type[b] == "RANK") {
      # RANK: all pairs observed, same set for every person
      items.b <- block.items[[b]]
      pairs.b <- t(combn(items.b, 2))
      for (p in 1:N.person) {
        for (pp in 1:n_pairs.b) {
          idx             <- idx + 1L
          Y[idx]          <- response[p, col_offset + pp]
          Idx_person[idx] <- p
          Idx_item_i[idx] <- pairs.b[pp, 1L]
          Idx_item_k[idx] <- pairs.b[pp, 2L]
          Idx_pair[idx]   <- match(paste0(pairs.b[pp, 1L], "_", pairs.b[pp, 2L]), pair.keys)
        }
      }
    } else {
      # MOLE / PICK: person-specific pairs
      for (p in 1:N.person) {
        pairs.cur <- pairs.value[[p]][[b]]
        if (!is.matrix(pairs.cur)) pairs.cur <- matrix(pairs.cur, ncol = 2)
        for (pp in 1:n_pairs.b) {
          idx             <- idx + 1L
          Y[idx]          <- response[p, col_offset + pp]
          Idx_person[idx] <- p
          Idx_item_i[idx] <- pairs.cur[pp, 1L]
          Idx_item_k[idx] <- pairs.cur[pp, 2L]
          Idx_pair[idx]   <- match(paste0(pairs.cur[pp, 1L], "_", pairs.cur[pp, 2L]), pair.keys)
        }
      }
    }
    col_offset <- col_offset + n_pairs.b
  }

  stopifnot(idx == N.obs, max(Idx_pair) <= N.pairs.all, min(Idx_pair) >= 1L)

  Idx_gamma_est <- sort(unique(Idx_pair))
  Idx_gamma_fix <- setdiff(seq_len(N.pairs.all), Idx_gamma_est)

  is_case_b <- (D == 2 && I.block == 2)
  if (is_case_b) {
    first_pair_items     <- block.items[[1]]
    Idx_lambda_fix       <- first_pair_items
    lambda_fix_val       <- rep(0.80, length(first_pair_items))

    Idx_lambda_est       <- setdiff(Idx_lambda_est_raw, Idx_lambda_fix)
  } else {
    Idx_lambda_est       <- Idx_lambda_est_raw
    Idx_lambda_fix       <- integer(0)
    lambda_fix_val       <- integer(0)
  }

  last_item_iblock <- unname(sapply(block.items, function(x) x[length(x)]))
  is_case_a        <- (I.block == 2 && D > 2)

  if (is_case_a || is_case_b) {
    Idx_psi_fix   <- 1:I.states
    psi_fix_val   <- rep(sqrt(0.5), I.states)
    Idx_psi_est   <- integer(0)
  } else {
    Idx_psi_fix   <- last_item_iblock
    psi_fix_val   <- rep(1.0, length(last_item_iblock))
    Idx_psi_est   <- setdiff(1:I.states, Idx_psi_fix)
  }

  Idx_psi_equal <- integer(0)
  Idx_psi_orig  <- integer(0)

  corr_df <- free_corr_npar(D)
  npar <- length(Idx_psi_est) + length(Idx_lambda_est) +
    length(Idx_gamma_est) + corr_df

  gamma_mu <- get_ctrl("gamma.mu", 0.0, control)
  gamma_fix_val <- rep(gamma_mu, length(Idx_gamma_fix))

  if (length(Idx_gamma_fix) > 0) {
    warning(sprintf(
      "TIRT: %d of %d possible pairwise comparisons were never observed (fc.type includes PICK/MOLE). The corresponding gamma parameters are fixed to the prior mean (%.2f) to avoid convergence issues from unobserved comparisons.",
      length(Idx_gamma_fix), N.pairs.all, gamma_mu
    ))
  }

  stan.data <- list(
    N_obs          = N.obs,
    N_pair         = N.pairs.all,
    N_person       = N.person,
    I_states       = I.states,
    D              = D,
    response       = Y,
    Idx_person     = Idx_person,
    Idx_item_i     = Idx_item_i,
    Idx_item_k     = Idx_item_k,
    Idx_pair       = Idx_pair,
    Q_matrix       = Q.matrix,
    N_lambda_est   = length(Idx_lambda_est),
    Idx_lambda_est = Idx_lambda_est,
    N_lambda_fix   = length(Idx_lambda_fix),
    Idx_lambda_fix = Idx_lambda_fix,
    lambda_fix_val = lambda_fix_val,
    N_psi_fix      = length(Idx_psi_fix),
    Idx_psi_fix    = Idx_psi_fix,
    psi_fix_val    = psi_fix_val,
    N_psi_est      = length(Idx_psi_est),
    Idx_psi_est    = Idx_psi_est,
    N_psi_equal    = length(Idx_psi_equal),
    Idx_psi_equal  = Idx_psi_equal,
    Idx_psi_orig   = Idx_psi_orig,
    N_gamma_fix    = length(Idx_gamma_fix),
    Idx_gamma_fix  = Idx_gamma_fix,
    gamma_fix_val  = gamma_fix_val,
    N_gamma_est    = length(Idx_gamma_est),
    Idx_gamma_est  = Idx_gamma_est,
    lambda_alpha   = get_ctrl("lambda.alpha", 0.40, control),
    lambda_beta    = get_ctrl("lambda.beta", 0.90, control),
    psi_mu         = get_ctrl("psi.mu", 0.80, control),
    psi_sigma      = get_ctrl("psi.sigma", 0.30, control),
    gamma_mu       = gamma_mu,
    gamma_sigma    = get_ctrl("gamma.sigma", 0.8, control),
    theta_mu       = get_ctrl("theta.mu", rep(0, D), control)
  )

  set.seed(common.method$seed)
  stan.obj <- suppressWarnings(
    rstan::sampling(
      object        = sm,
      data          = stan.data,
      chains        = chains,
      iter          = iter,
      warmup        = warmup,
      thin          = thin,
      init          = init,
      algorithm     = algorithm,
      cores         = mcmc_env$cores,
      control       = stan_control,
      verbose       = common.method$vis,
      show_messages = common.method$vis,
      open_progress = FALSE
    )
  )

  MCMC.obj <- rstan::extract(stan.obj, permuted = TRUE)
  stan.sum <- rstan::summary(stan.obj)$summary

  theta_est <- extract_theta_stan(MCMC.obj, stan.sum, N.person, D, rownames(response))

  lambda_raw  <- apply(MCMC.obj$lambda, 2, mean)
  lambda.se   <- apply(MCMC.obj$lambda, 2, sd)
  lambda.Rhat <- rep(NA_real_, I.states)
  lam.idx     <- grep("^lambda\\[", rownames(stan.sum))
  if (length(lam.idx) > 0) lambda.Rhat <- stan.sum[lam.idx, "Rhat"]

  lambda.free <- rep(FALSE, I.states)
  lambda.free[Idx_lambda_est] <- TRUE
  lambda.se[!lambda.free] <- NA_real_
  lambda.Rhat[!lambda.free] <- NA_real_

  lambda_sign <- apply(Q.matrix, 1, function(x) {
    nz <- which(x != 0)
    if (length(nz) == 0) return(1)
    sign(x[nz[1]])
  })
  lambda <- lambda_raw * lambda_sign

  psi_sd      <- apply(MCMC.obj$psi, 2, mean)
  psi_sd.se   <- apply(MCMC.obj$psi, 2, sd)
  psi2        <- psi_sd^2
  psi2.se     <- 2 * psi_sd * psi_sd.se
  psi2.Rhat   <- rep(NA_real_, I.states)
  psi.idx     <- grep("^psi\\[", rownames(stan.sum))
  if (length(psi.idx) > 0) psi2.Rhat <- stan.sum[psi.idx, "Rhat"]

  psi.free <- rep(FALSE, I.states)
  psi.free[Idx_psi_est] <- TRUE
  psi_sd.se[!psi.free] <- NA_real_
  psi2.se[!psi.free] <- NA_real_
  psi2.Rhat[!psi.free] <- NA_real_

  gamma      <- apply(MCMC.obj$gamma, 2, mean)
  gamma.se   <- apply(MCMC.obj$gamma, 2, sd)
  gamma.Rhat <- rep(NA_real_, N.pairs.all)
  gam.idx    <- grep("^gamma\\[", rownames(stan.sum))
  if (length(gam.idx) > 0) gamma.Rhat <- stan.sum[gam.idx, "Rhat"]

  gamma.free <- rep(FALSE, N.pairs.all)
  gamma.free[Idx_gamma_est] <- TRUE
  gamma.se[!gamma.free] <- NA_real_
  gamma.Rhat[!gamma.free] <- NA_real_

  gamma.matrix      <- matrix(NA_real_, nrow = I.states, ncol = I.states)
  gamma.matrix.se   <- matrix(NA_real_, nrow = I.states, ncol = I.states)
  gamma.matrix.Rhat <- matrix(NA_real_, nrow = I.states, ncol = I.states)
  for (r in seq_len(N.pairs.all)) {
    i <- pairs.matrix[r, 1L]; k <- pairs.matrix[r, 2L]
    gamma.matrix[i, k]      <- gamma[r]
    gamma.matrix.se[i, k]   <- gamma.se[r]
    gamma.matrix.Rhat[i, k] <- gamma.Rhat[r]
    gamma.matrix[k, i]      <- -gamma[r]
    gamma.matrix.se[k, i]   <- gamma.se[r]
    gamma.matrix.Rhat[k, i] <- gamma.Rhat[r]
  }

  Corr_est <- extract_corr_stan(MCMC.obj, stan.sum, D)
  par_est <- list(
    est   = cbind(lambda * abs(Q.matrix), psi2),
    se    = cbind(lambda.se * abs(Q.matrix), psi2.se),
    Rhat  = cbind(lambda.Rhat * abs(Q.matrix), psi2.Rhat),
    free  = tirt_par_free_mask(Q.matrix, lambda.free, psi.free)
  )

  colnames(par_est$est) <- colnames(par_est$se) <- colnames(par_est$Rhat) <-
    colnames(par_est$free) <- c(paste0("Dim.", 1:D), "psi2")
  rownames(par_est$est) <- rownames(par_est$se) <- rownames(par_est$Rhat) <-
    rownames(par_est$free) <- paste0("item", seq_len(I.states))
  par_est$se[!par_est$free] <- NA_real_
  par_est$Rhat[!par_est$free] <- NA_real_

  gamma_est <- list(est=gamma, se=gamma.se, Rhat=gamma.Rhat, free=gamma.free)
  gamma.matrix_est <- list(est = gamma.matrix, se = gamma.matrix.se, Rhat = gamma.matrix.Rhat)

  results <- list(
    npar         = npar,
    method       = "stan",
    theta        = theta_est,
    par          = par_est,
    gamma        = gamma_est,
    Corr         = Corr_est,
    gamma.matrix = gamma.matrix_est,
    Q.matrix     = Q.matrix,
    block.items  = block.items,
    pairs.matrix = pairs.matrix,
    pairs.value  = pairs.value,
    response     = response,
    fc.type      = fc.type,
    stan.obj     = stan.obj,
    MCMC.obj     = MCMC.obj,
    log_lik      = MCMC.obj$log_lik,
    call = call,
    arguments = list(
      data        = data,
      D           = D,
      Q.matrix    = Q.matrix,
      method      = "stan",
      block.items = block.items,
      fc.type     = fc.type,
      chains      = chains,
      iter        = iter,
      warmup      = warmup,
      thin        = thin,
      init        = init,
      algorithm   = algorithm,
      cores       = common.method$cores,
      vis         = common.method$vis,
      seed        = common.method$seed,
      control.model = control,
      control.method = fit_effective_method_control(
        control.method, method = stan.method, common = common.method
      )
    )
  )
  class(results) <- "TIRT"

  L <- get_ctrl("L", NULL, control)
  results$logLik <- logLik(results, L=L)

  return(results)
}
