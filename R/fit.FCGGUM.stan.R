#' Fit FCGGUM Model via Stan
#'
#' Fits the Forced-Choice Generalized Graded Unfolding Model (FCGGUM) using
#' Bayesian MCMC via Stan.  The FCGGUM is a forced-choice wrapper around the
#' MGGUM: item-level endorsement probabilities are obtained from the MGGUM
#' binary endorsement probability, transformed to logits, and full rankings
#' are modeled as sequential Luce/PICK decisions over the remaining
#' statements in each block.
#'
#' @param data An N x B data frame where each column is a block and each cell
#'   contains the ranking order as a character string (e.g., "3>1>2"), or
#'   alternatively an N x B integer matrix of pattern indices (1-based) as
#'   returned by \code{sim.data.FCGGUM()$response}.
#' @param Q.matrix An I x D matrix with values -1, 0, or 1.
#'   \itemize{
#'     \item \code{1}: active dimension with \code{delta} constrained/regularized on the positive side
#'     \item \code{-1}: active dimension with \code{delta} constrained/regularized on the negative side
#'     \item \code{0}: \code{a} and \code{delta} fixed to 0 (inactive dimension)
#'   }
#'   The sign controls the side of the item-location prior/support, not the
#'   sign of discrimination; \code{a} is nonnegative whenever \code{|Q| = 1}.
#'   If \code{NULL}, a default matrix of all 1's is used.
#' @param block.items A list of length B; each element is an integer vector of
#'   (global) item indices belonging to that block.  If \code{NULL}, the block
#'   structure is auto-detected only from full-rank \code{"RANK"} string
#'   \code{data}; it is required when \code{data} is an integer matrix of
#'   pattern indices or when \code{fc.type} contains \code{"MOLE"} or
#'   \code{"PICK"}.
#' @param D Number of latent dimensions (must be >= 2).  If \code{NULL},
#'   \code{D} is derived from \code{ncol(Q.matrix)} when \code{Q.matrix}
#'   is provided.
#' @param control.model Optional named list of model-level controls and
#'   hyperparameters. Supported prior entries include \code{a.mu},
#'   \code{a.sigma}, \code{delta.pos.mu}, \code{delta.pos.sigma},
#'   \code{delta.neg.mu}, \code{delta.neg.sigma}, \code{tau.mu},
#'   \code{tau.sigma}, and \code{theta.mu}.
#' @param control.method Optional named list of method-specific controls. For
#'   Stan this includes \code{cores}, \code{vis}, \code{seed}, \code{chains},
#'   \code{iter}, \code{warmup}, \code{thin}, \code{init}, \code{algorithm},
#'   and Stan sampler controls.
#'
#' @noRd
fit.FCGGUM.stan <- function(data, Q.matrix = NULL, block.items = NULL, D = NULL, fc.type = "RANK",
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

  sm <- stanmodels[["FCGGUM"]]

  if (is.null(sm)) {
    stop("Stan model FCGGUM not found !", call. = FALSE)
  }

  mcmc_env <- setup_mcmc_env(common.method$cores)
  on.exit(restore_mcmc_env(mcmc_env), add = TRUE)

  stan_control <- build_stan_control(algorithm, control.method)

  # ---- Validate fc.type ----
  fc.type.valid <- c("RANK", "MOLE", "PICK")
  if (!all(fc.type %in% fc.type.valid)) {
    stop("'fc.type' must be one of: ", paste(fc.type.valid, collapse = ", "), call. = FALSE)
  }

  data.matrix <- as.matrix(data)
  N.block.data <- ncol(data.matrix)
  if (is.null(N.block.data) || N.block.data < 1L) {
    stop("'data' must have at least one column/block.", call. = FALSE)
  }
  if (length(fc.type) != N.block.data) {
    fc.type <- rep(fc.type[1], N.block.data)
  }

  # ---- Auto-detect block items from full-rank string data if not provided ----
  is.pattern.index.data <- is.numeric(data.matrix) || is.integer(data.matrix)

  if (is.null(block.items)) {
    if (is.pattern.index.data) {
      stop("'block.items' is required when 'data' is an integer matrix of pattern indices.", call. = FALSE)
    }
    if (any(fc.type != "RANK")) {
      stop("'block.items' is required for MOLE/PICK data because partial rankings do not identify all block items.", call. = FALSE)
    }
    block.items <- get.block.items.from.data(data)
  }

  # ---- Convert ranking string data to response (pattern indices) ----
  if (is.pattern.index.data) {
    response <- matrix(as.integer(data.matrix), nrow = nrow(data.matrix), ncol = ncol(data.matrix),
                       dimnames = dimnames(data.matrix))
  } else {
    response <- get.response.from.data(data, block.items, fc.type = fc.type)
  }

  # ---- Validate block.items ----
  if (!is.list(block.items)) {
    stop("'block.items' must be a list.", call. = FALSE)
  }
  N.block <- length(block.items)

  all_items <- unlist(block.items)
  if (anyDuplicated(all_items)) {
    stop("'block.items' contains duplicate items.", call. = FALSE)
  }
  I <- length(all_items)
  if (!setequal(all_items, seq_len(I))) {
    stop("'block.items' must contain exactly the item indices 1 through ", I, ".", call. = FALSE)
  }
  if (any(vapply(block.items, length, integer(1L)) < 2L)) {
    stop("Each block in 'block.items' must contain at least two items.", call. = FALSE)
  }

  # ---- Validate response ----
  response <- as.matrix(response)
  N <- nrow(response)
  if (ncol(response) != N.block) {
    stop("'response' must have ", N.block, " columns (one per block).", call. = FALSE)
  }

  if (length(fc.type) != N.block) {
    fc.type <- rep(fc.type[1], N.block)
  }

  # ---- Map fc.type to integer codes for Stan ----
  fc.type.int <- match(fc.type, fc.type.valid)

  # ---- Item categories (binary agree/disagree for FCGGUM) ----
  length.poly <- rep(2L, I)
  max_poly <- 2L

  # ---- Validate Q.matrix and derive D ----
  if (!is.null(Q.matrix)) {
    if (!is.matrix(Q.matrix) || anyNA(Q.matrix) || !all(Q.matrix %in% c(-1, 0, 1))) {
      stop("'Q.matrix' must be a matrix containing only -1, 0, or 1", call. = FALSE)
    }
    D.q <- ncol(Q.matrix)
    if (!is.null(D) && D != D.q) {
      warning("'D' (", D, ") does not match ncol(Q.matrix) (", D.q, "); using ncol(Q.matrix).")
    }
    D <- D.q
    if (nrow(Q.matrix) != I) {
      stop("'Q.matrix' dimensions (", nrow(Q.matrix), "x", ncol(Q.matrix),
           ") do not match data (", I, "x", D, ").", call. = FALSE)
    } else if (any(rowSums(abs(Q.matrix)) < 1)) {
      stop("Each item must measure at least one trait (each row of 'Q.matrix' needs a non-zero entry)", call. = FALSE)
    }
  }

  # ---- Validate D ----
  if (is.null(D) || is.na(D)) {
    stop("Cannot determine 'D'. Please provide either 'Q.matrix' or 'D'.", call. = FALSE)
  }
  if (D < 2) {
    stop("FCGGUM requires at least 2 latent dimensions (D >= 2).", call. = FALSE)
  }

  if (is.null(Q.matrix)) {
    Q.matrix <- matrix(1, I, D)
  }

  corr_df <- free_corr_npar(D)

  # ---- Precompute patterns for each block ----
  patterns.total <- vector("list", N.block)
  for (b in 1:N.block) {
    patterns.total[[b]] <- get_permutations(block.items[[b]],
                                            length(block.items[[b]]))
  }

  patterns.obs <- vector("list", N.block)
  for (b in 1:N.block) {
    patterns.obs[[b]] <- get_permutations(
      block.items[[b]],
      switch(fc.type[b],
             RANK = length(block.items[[b]]),
             MOLE = 2,
             PICK = 1)
    )
  }

  for (b in 1:N.block) {
    if (anyNA(response[, b]) ||
        any(response[, b] < 1L | response[, b] > nrow(patterns.obs[[b]]))) {
      stop("'data' contains invalid pattern indices for block ", b,
           "; valid values are 1 through ", nrow(patterns.obs[[b]]), ".", call. = FALSE)
    }
  }

  # ---- Build flattened arrays for Stan ----
  n_items     <- integer(N.block)
  n_total     <- integer(N.block)
  n_obs       <- integer(N.block)
  fc_len      <- integer(N.block)

  block_items_flat   <- integer(0)
  patterns_total_flat <- integer(0)
  patterns_obs_flat   <- integer(0)

  item_start  <- integer(N.block + 1L)
  total_start <- integer(N.block + 1L)
  obs_start   <- integer(N.block + 1L)

  item_start[1]  <- 1L
  total_start[1] <- 1L
  obs_start[1]   <- 1L

  for (b in 1:N.block) {
    items.b <- block.items[[b]]
    K <- length(items.b)
    n_items[b] <- K

    fc_len[b] <- switch(fc.type[b],
                        RANK = K,
                        MOLE = 2L,
                        PICK = 1L)

    # Full ranking patterns to LOCAL indices 1..K
    pat.total <- patterns.total[[b]]
    n_total[b] <- nrow(pat.total)
    pat.total.local <- matrix(match(pat.total, items.b), nrow = n_total[b])
    patterns_total_flat <- c(patterns_total_flat, as.integer(t(pat.total.local)))

    # FC-type patterns to LOCAL indices 1..K
    pat.obs <- patterns.obs[[b]]
    n_obs[b] <- nrow(pat.obs)
    pat.obs.local <- matrix(match(pat.obs, items.b), nrow = n_obs[b])
    patterns_obs_flat <- c(patterns_obs_flat, as.integer(t(pat.obs.local)))

    # Global block items
    block_items_flat <- c(block_items_flat, items.b)

    item_start[b + 1]  <- item_start[b] + K
    total_start[b + 1] <- total_start[b] + n_total[b] * K
    obs_start[b + 1]   <- obs_start[b] + n_obs[b] * fc_len[b]
  }

  total_block_items <- length(block_items_flat)
  total_cells_total <- length(patterns_total_flat)
  total_cells_obs   <- length(patterns_obs_flat)

  # ---- Convert response to long format ----
  N_obs    <- N * N.block
  pid_long <- rep(1:N, times = N.block)
  bid_long <- rep(1:N.block, each = N)
  y_long   <- as.integer(response)

  # ---- Stan data ----
  theta_mu_raw <- get_ctrl("theta.mu", rep(0, D), control)
  theta_mu     <- array(theta_mu_raw, dim = D)

  stan.data <- list(
    N_obs              = N_obs,
    N                  = N,
    B                  = N.block,
    I                  = I,
    D                  = D,
    max_poly           = max_poly,
    pid                = pid_long,
    bid                = bid_long,
    y                  = y_long,
    length_poly        = as.integer(length.poly),
    Q_matrix           = Q.matrix,
    n_items            = n_items,
    fc_type            = fc.type.int,
    total_block_items  = total_block_items,
    block_items        = block_items_flat,
    item_start         = item_start,
    n_total            = n_total,
    total_cells_total  = total_cells_total,
    patterns_total     = patterns_total_flat,
    total_start        = total_start,
    n_obs              = n_obs,
    total_cells_obs    = total_cells_obs,
    patterns_obs       = patterns_obs_flat,
    obs_start          = obs_start,
    fc_len             = fc_len,
    a_mu               = get_ctrl("a.mu", 0, control),
    a_sigma            = get_ctrl("a.sigma", 0.5, control),
    delta_pos_mu       = get_ctrl("delta.pos.mu", 1, control),
    delta_pos_sigma    = get_ctrl("delta.pos.sigma", 0.5, control),
    delta_neg_mu       = get_ctrl("delta.neg.mu", -1, control),
    delta_neg_sigma    = get_ctrl("delta.neg.sigma", 0.5, control),
    tau_mu             = get_ctrl("tau.mu", 0, control),
    tau_sigma          = get_ctrl("tau.sigma", 0.5, control),
    theta_mu           = theta_mu
  )

  set.seed(common.method$seed)

  # ---- Run Stan ----
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
      show_messages = common.method$vis
    )
  )

  MCMC.obj <- rstan::extract(stan.obj, permuted = TRUE)
  stan.sum <- rstan::summary(stan.obj)$summary

  # ---- Extract theta ----
  theta_est <- extract_theta_stan(MCMC.obj, stan.sum, N, D, rownames(response))

  # ---- Build par matrix: I x (D + D + max_poly) ----
  par      <- matrix(NA_real_, nrow = I, ncol = D + D + max_poly)
  par.se   <- matrix(NA_real_, nrow = I, ncol = D + D + max_poly)
  par.Rhat <- matrix(NA_real_, nrow = I, ncol = D + D + max_poly)

  a_arr <- ensure_3d(MCMC.obj$a)
  a     <- apply(a_arr, c(2, 3), mean)
  a.se  <- apply(a_arr, c(2, 3), sd)

  par[, 1:D]      <- a
  par.se[, 1:D]   <- a.se
  par.Rhat[, 1:D] <- extract_rhat_matrix(stan.sum, "a", I, D)

  delta_arr <- ensure_3d(MCMC.obj$delta)
  delta     <- apply(delta_arr, c(2, 3), mean)
  delta.se  <- apply(delta_arr, c(2, 3), sd)

  par[, (D + 1):(D + D)]      <- delta
  par.se[, (D + 1):(D + D)]   <- delta.se
  par.Rhat[, (D + 1):(D + D)] <- extract_rhat_matrix(stan.sum, "delta", I, D)

  par[, D + D + 1]      <- 0
  par.se[, D + D + 1]   <- NA_real_
  par.Rhat[, D + D + 1] <- NA_real_

  tau_arr <- ensure_3d(MCMC.obj$tau)
  tau     <- apply(tau_arr, c(2, 3), mean)
  tau.se  <- apply(tau_arr, c(2, 3), sd)
  for (i in 1:I) {
    Ki <- length.poly[i]
    if (Ki > 1) {
      tau_cols <- (D + D + 2):(D + D + Ki)
      par[i, tau_cols]      <- tau[i, 2:Ki]
      par.se[i, tau_cols]   <- tau.se[i, 2:Ki]
      tau_idx_i <- grep(paste0("^tau\\[", i, ","), rownames(stan.sum))
      tau_rhat_i <- stan.sum[tau_idx_i, "Rhat"]
      par.Rhat[i, tau_cols] <- tau_rhat_i[2:Ki]
    }
  }

  item_names <- paste0("item", all_items)
  par.free <- mggum_par_free_mask(Q.matrix, length.poly, max_poly)
  colnames(par) <- colnames(par.se) <- colnames(par.Rhat) <- colnames(par.free) <-
    c(paste0("a", 1:D),
      paste0("delta", 1:D),
      paste0("tau", 0:(max_poly - 1)))
  rownames(par) <- rownames(par.se) <- rownames(par.Rhat) <- rownames(par.free) <-
    item_names
  par.se[!par.free] <- NA_real_
  par.Rhat[!par.free] <- NA_real_
  npar <- sum(par.free) + corr_df

  # ---- Extract Corr ----
  Corr_est <- extract_corr_stan(MCMC.obj, stan.sum, D)

  # ---- Assemble results ----
  par_est   <- list(est = par, se = par.se, Rhat = par.Rhat, free = par.free)

  results <- list(
    npar           = npar,
    method         = "stan",
    theta          = theta_est,
    par            = par_est,
    Corr           = Corr_est,
    stan.obj       = stan.obj,
    MCMC.obj       = MCMC.obj,
    log_lik         = MCMC.obj$log_lik,
    Q.matrix       = Q.matrix,
    block.items    = block.items,
    fc.type        = fc.type,
    response       = response,
    length.poly    = length.poly,
    patterns       = patterns.obs,
    patterns.total = patterns.total,
    call           = call,
    arguments      = list(
      data        = data,
      Q.matrix    = Q.matrix,
      method      = "stan",
      block.items = block.items,
      D           = D,
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

  class(results) <- "FCGGUM"

  L <- get_ctrl("L", NULL, control)
  results$logLik <- logLik.FCGGUM(results, L = L)

  return(results)
}
