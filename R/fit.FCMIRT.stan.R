#' Fit FCMIRT Model via Stan
#'
#' Fits the forced-choice MIRT model using Bayesian MCMC via Stan. Item-level
#' endorsement probabilities are defined by an MIRT model
#' (\code{"m1pl"}, \code{"m2pl"}, \code{"m3pl"}, or \code{"m4pl"}) and block
#' responses are modeled through the same RANK, MOLE, and PICK mechanisms used
#' by \code{\link[ForceChoice]{sim.data.FCMIRT}}.
#' The block difficulty/intercept parameters \code{b} are identified by a
#' zero-sum constraint within each forced-choice block; Stan samples the
#' \code{K - 1} free difficulties and reconstructs the final one.
#'
#' @param data An N x B data frame where each column is a forced-choice block
#'   and each cell is a ranking string such as \code{"3>1>2"}, or an N x B
#'   integer matrix of 1-based pattern indices.
#' @param model MIRT model type: \code{"m1pl"}, \code{"m2pl"}, \code{"m3pl"},
#'   or \code{"m4pl"}.
#' @param Q.matrix An I x D 0/1 matrix indicating active item loadings. For
#'   \code{"m1pl"}, active slopes are fixed to 1 and inactive slopes to 0.
#' @param block.items A list of length B giving global item indices in each
#'   block. Required for integer pattern-index data and for MOLE/PICK data.
#' @param D Number of latent dimensions. If \code{Q.matrix} is supplied,
#'   \code{D} is reset to \code{ncol(Q.matrix)}.
#' @param fc.type Forced-choice response type: \code{"RANK"}, \code{"MOLE"}, or
#'   \code{"PICK"}. A scalar value is recycled to all blocks.
#' @param control.model Optional named list of model-level controls and
#'   hyperparameters. Supported prior entries include \code{a.mu},
#'   \code{a.sigma}, \code{b.mu}, \code{b.sigma}, \code{c.mu},
#'   \code{c.sigma}, \code{d.mu}, \code{d.sigma}, and \code{theta.mu}. The
#'   Stan default for \code{a.sigma} is deliberately weak because
#'   forced-choice blocks identify relative utilities and tight common slope
#'   priors can over-shrink item slopes within the same block.
#' @param control.method Optional named list of method-specific controls. For
#'   Stan this includes \code{cores}, \code{vis}, \code{seed}, \code{chains},
#'   \code{iter}, \code{warmup}, \code{thin}, \code{init}, \code{algorithm},
#'   and Stan sampler controls.
#'
#' @noRd
fit.FCMIRT.stan <- function(data, model = "m2pl", Q.matrix = NULL, block.items = NULL,
                       D = NULL, fc.type = "RANK",
                       control.model = NULL,
                       control.method = NULL,
                       .call = NULL) {
  call <- if (is.null(.call)) match.call() else .call
  control <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  stan.method <- fc_stan_method_control(control.method)

  model <- tolower(model)
  model_type <- switch(model,
    "m1pl" = 1L,
    "m2pl" = 2L,
    "m3pl" = 3L,
    "m4pl" = 4L,
    stop("Unknown FCMIRT model type: ", model, call. = FALSE)
  )
  chains <- stan.method$chains
  iter <- stan.method$iter
  warmup <- stan.method$warmup
  thin <- stan.method$thin
  init <- stan.method$init
  algorithm <- stan.method$algorithm

  sm <- stanmodels[["FCMIRT"]]
  if (is.null(sm)) {
    stop("Stan model FCMIRT not found !", call. = FALSE)
  }

  mcmc_env <- setup_mcmc_env(common.method$cores)
  on.exit(restore_mcmc_env(mcmc_env), add = TRUE)

  stan_control <- build_stan_control(algorithm, control.method)

  fc.type.valid <- c("RANK", "MOLE", "PICK")
  fc.type <- toupper(fc.type)
  if (!all(fc.type %in% fc.type.valid)) {
    stop("'fc.type' must be one of: ", paste(fc.type.valid, collapse = ", "), call. = FALSE)
  }

  data.matrix <- as.matrix(data)
  N.block.data <- ncol(data.matrix)
  if (is.null(N.block.data) || N.block.data < 1L) {
    stop("'data' must have at least one column/block.", call. = FALSE)
  }
  if (length(fc.type) != N.block.data) {
    fc.type <- rep(fc.type[1L], N.block.data)
  }

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

  if (is.pattern.index.data) {
    response <- matrix(as.integer(data.matrix),
                       nrow = nrow(data.matrix), ncol = ncol(data.matrix),
                       dimnames = dimnames(data.matrix))
  } else {
    response <- get.response.from.data(data, block.items, fc.type = fc.type)
  }

  if (!is.list(block.items)) {
    stop("'block.items' must be a list.", call. = FALSE)
  }
  N.block <- length(block.items)

  all_items <- unlist(block.items, use.names = FALSE)
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

  response <- as.matrix(response)
  N <- nrow(response)
  if (ncol(response) != N.block) {
    stop("'response' must have ", N.block, " columns (one per block).", call. = FALSE)
  }

  if (length(fc.type) != N.block) {
    fc.type <- rep(fc.type[1L], N.block)
  }
  fc.type.int <- match(fc.type, fc.type.valid)

  if (!is.null(Q.matrix)) {
    Q.matrix <- as.matrix(Q.matrix)
    if (anyNA(Q.matrix) || !all(Q.matrix %in% c(0, 1))) {
      stop("'Q.matrix' must be a matrix containing only 0 or 1.", call. = FALSE)
    }
    if (nrow(Q.matrix) != I) {
      stop("'Q.matrix' must have ", I, " rows, one for each item.", call. = FALSE)
    }
    if (any(rowSums(Q.matrix) < 1L)) {
      stop("Each item must measure at least one trait (each row of 'Q.matrix' needs a 1).", call. = FALSE)
    }
    D.q <- ncol(Q.matrix)
    if (!is.null(D) && D != D.q) {
      warning("'D' (", D, ") does not match ncol(Q.matrix) (", D.q, "); using ncol(Q.matrix).")
    }
    D <- D.q
    storage.mode(Q.matrix) <- "numeric"
  }

  if (is.null(D) || is.na(D)) {
    stop("Cannot determine 'D'. Please provide either 'Q.matrix' or 'D'.", call. = FALSE)
  }
  D <- as.integer(D)
  if (D < 1L) {
    stop("'D' must be a positive integer.", call. = FALSE)
  }

  if (is.null(Q.matrix)) {
    Q.matrix <- matrix(1, I, D)
  }

  corr_df <- free_corr_npar(D)

  patterns.total <- vector("list", N.block)
  patterns.obs <- vector("list", N.block)
  for (b in seq_len(N.block)) {
    patterns.total[[b]] <- get_permutations(block.items[[b]],
                                            length(block.items[[b]]))
    patterns.obs[[b]] <- get_permutations(
      block.items[[b]],
      switch(fc.type[b],
             RANK = length(block.items[[b]]),
             MOLE = 2L,
             PICK = 1L)
    )
  }

  for (b in seq_len(N.block)) {
    if (anyNA(response[, b]) ||
        any(response[, b] < 1L | response[, b] > nrow(patterns.obs[[b]]))) {
      stop("'data' contains invalid pattern indices for block ", b,
           "; valid values are 1 through ", nrow(patterns.obs[[b]]), ".", call. = FALSE)
    }
  }

  n_items <- integer(N.block)
  n_total <- integer(N.block)
  n_obs <- integer(N.block)
  fc_len <- integer(N.block)

  block_items_flat <- integer(0)
  patterns_total_flat <- integer(0)
  patterns_obs_flat <- integer(0)

  item_start <- integer(N.block + 1L)
  total_start <- integer(N.block + 1L)
  obs_start <- integer(N.block + 1L)
  item_start[1L] <- 1L
  total_start[1L] <- 1L
  obs_start[1L] <- 1L

  for (b in seq_len(N.block)) {
    items.b <- block.items[[b]]
    K <- length(items.b)
    n_items[b] <- K
    fc_len[b] <- switch(fc.type[b], RANK = K, MOLE = 2L, PICK = 1L)

    pat.total <- patterns.total[[b]]
    n_total[b] <- nrow(pat.total)
    pat.total.local <- matrix(match(pat.total, items.b), nrow = n_total[b])
    patterns_total_flat <- c(patterns_total_flat, as.integer(t(pat.total.local)))

    pat.obs <- patterns.obs[[b]]
    n_obs[b] <- nrow(pat.obs)
    pat.obs.local <- matrix(match(pat.obs, items.b), nrow = n_obs[b])
    patterns_obs_flat <- c(patterns_obs_flat, as.integer(t(pat.obs.local)))

    block_items_flat <- c(block_items_flat, items.b)

    item_start[b + 1L] <- item_start[b] + K
    total_start[b + 1L] <- total_start[b] + n_total[b] * K
    obs_start[b + 1L] <- obs_start[b] + n_obs[b] * fc_len[b]
  }

  N_obs <- N * N.block
  pid_long <- rep(seq_len(N), times = N.block)
  bid_long <- rep(seq_len(N.block), each = N)
  y_long <- as.integer(response)

  theta_mu_raw <- get_ctrl("theta.mu", rep(0, D), control)
  theta_mu <- array(theta_mu_raw, dim = D)

  stan.data <- list(
    N_obs = N_obs,
    N = N,
    B = N.block,
    I = I,
    D = D,
    model_type = model_type,
    pid = pid_long,
    bid = bid_long,
    y = y_long,
    Q_matrix = Q.matrix,
    n_items = n_items,
    fc_type = fc.type.int,
    total_block_items = length(block_items_flat),
    block_items = block_items_flat,
    item_start = item_start,
    n_total = n_total,
    total_cells_total = length(patterns_total_flat),
    patterns_total = patterns_total_flat,
    total_start = total_start,
    n_obs = n_obs,
    total_cells_obs = length(patterns_obs_flat),
    patterns_obs = patterns_obs_flat,
    obs_start = obs_start,
    fc_len = fc_len,
    a_mu = get_ctrl("a.mu", 0.0, control),
    a_sigma = get_ctrl("a.sigma", 1.0, control),
    b_mu = get_ctrl("b.mu", 0.0, control),
    b_sigma = get_ctrl("b.sigma", 1.0, control),
    c_mu = get_ctrl("c.mu", 0.0, control),
    c_sigma = get_ctrl("c.sigma", 0.35, control),
    d_mu = get_ctrl("d.mu", 0.65, control),
    d_sigma = get_ctrl("d.sigma", 1.0, control),
    theta_mu = theta_mu
  )

  set.seed(common.method$seed)

  stan.obj <- suppressWarnings(
    rstan::sampling(
      object = sm,
      data = stan.data,
      chains = chains,
      iter = iter,
      warmup = warmup,
      thin = thin,
      init = init,
      algorithm = algorithm,
      cores = mcmc_env$cores,
      control = stan_control,
      verbose = common.method$vis,
      show_messages = common.method$vis
    )
  )

  MCMC.obj <- rstan::extract(stan.obj, permuted = TRUE)
  stan.sum <- rstan::summary(stan.obj)$summary

  theta_est <- extract_theta_stan(MCMC.obj, stan.sum, N, D, rownames(response))

  par <- matrix(NA_real_, nrow = I, ncol = D + 3L)
  par.se <- matrix(NA_real_, nrow = I, ncol = D + 3L)
  par.Rhat <- matrix(NA_real_, nrow = I, ncol = D + 3L)
  par.free <- fcmirt_par_free_mask(model_type, Q.matrix, block.items)

  if (model == "m1pl") {
    par[, seq_len(D)] <- Q.matrix
    par.se[, seq_len(D)] <- NA_real_
    par.Rhat[, seq_len(D)] <- NA_real_
  } else {
    a_arr <- ensure_3d(MCMC.obj$a)
    a <- apply(a_arr, c(2, 3), mean)
    a.se <- apply(a_arr, c(2, 3), sd)
    par[, seq_len(D)] <- a
    par.se[, seq_len(D)] <- a.se
    par.Rhat[, seq_len(D)] <- extract_rhat_matrix(stan.sum, "a", I, D)
  }

  b_est <- apply(MCMC.obj$b, 2, mean)
  par[, D + 1L] <- b_est
  b_free_items <- unlist(lapply(block.items, function(items) {
    items[-length(items)]
  }), use.names = FALSE)
  b_free_draws <- as.matrix(MCMC.obj$b_free)
  par.se[b_free_items, D + 1L] <- apply(b_free_draws, 2L, sd)
  par.Rhat[b_free_items, D + 1L] <-
    extract_rhat_vector(stan.sum, "b_free", length(b_free_items))

  if (model %in% c("m3pl", "m4pl")) {
    c_est <- apply(MCMC.obj$c, 2, mean)
    c.se <- apply(MCMC.obj$c, 2, sd)
    c.idx <- grep("^c\\[", rownames(stan.sum))
    par[, D + 2L] <- c_est
    par.se[, D + 2L] <- c.se
    par.Rhat[, D + 2L] <- stan.sum[c.idx, "Rhat"]
  } else {
    par[, D + 2L] <- 0
    par.se[, D + 2L] <- NA_real_
    par.Rhat[, D + 2L] <- NA_real_
  }

  if (model == "m4pl") {
    d_est <- apply(MCMC.obj$d, 2, mean)
    d.se <- apply(MCMC.obj$d, 2, sd)
    d.idx <- grep("^d\\[", rownames(stan.sum))
    par[, D + 3L] <- d_est
    par.se[, D + 3L] <- d.se
    par.Rhat[, D + 3L] <- stan.sum[d.idx, "Rhat"]
  } else {
    par[, D + 3L] <- 1
    par.se[, D + 3L] <- NA_real_
    par.Rhat[, D + 3L] <- NA_real_
  }

  item_names <- paste0("item", all_items)
  colnames(par) <- colnames(par.se) <- colnames(par.Rhat) <- colnames(par.free) <-
    c(paste0("a", seq_len(D)), "b", "c", "d")
  rownames(par) <- rownames(par.se) <- rownames(par.Rhat) <- rownames(par.free) <-
    item_names
  par.se[!par.free] <- NA_real_
  par.Rhat[!par.free] <- NA_real_
  npar <- sum(par.free) + corr_df

  Corr_est <- extract_corr_stan(MCMC.obj, stan.sum, D)

  results <- list(
    npar = npar,
    method = "stan",
    theta = theta_est,
    par = list(est = par, se = par.se, Rhat = par.Rhat, free = par.free),
    Corr = Corr_est,
    stan.obj = stan.obj,
    MCMC.obj = MCMC.obj,
    log_lik = MCMC.obj$log_lik,
    model = model,
    Q.matrix = Q.matrix,
    block.items = block.items,
    fc.type = fc.type,
    response = response,
    patterns = patterns.obs,
    patterns.total = patterns.total,
    call = call,
    arguments = list(
      data = data,
      model = model,
      Q.matrix = Q.matrix,
      method = "stan",
      block.items = block.items,
      D = D,
      fc.type = fc.type,
      chains = chains,
      iter = iter,
      warmup = warmup,
      thin = thin,
      init = init,
      algorithm = algorithm,
      cores = common.method$cores,
      vis = common.method$vis,
      seed = common.method$seed,
      control.model = control,
      control.method = fit_effective_method_control(
        control.method, method = stan.method, common = common.method
      )
    )
  )

  class(results) <- "FCMIRT"

  L <- get_ctrl("L", NULL, control)
  results$logLik <- logLik.FCMIRT(results, L = L)

  return(results)
}
