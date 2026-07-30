#' Fit MGPCM Model via Stan
#'
#' @param response An N x I matrix of integer responses (0, 1, ..., K_i - 1)
#' @param D Number of latent dimensions
#' @param Q.matrix An I x D 0/1 matrix indicating which discrimination parameters
#'   are freely estimated (1) and which are fixed to 0 (default: NULL, builds
#'   a triangular identification structure).
#' @param length.poly Optional category counts for the items.
#' @param control.model Optional named list of model-level controls and
#'   hyperparameters.
#' @param control.method Optional named list of method-specific controls. For
#'   Stan this includes \code{cores}, \code{vis}, \code{seed}, \code{chains},
#'   \code{iter}, \code{warmup}, \code{thin}, \code{init}, \code{algorithm},
#'   and Stan sampler controls.
#'
#' @noRd
fit.MGPCM.stan <- function(response, D = 2, Q.matrix = NULL,
                       length.poly = NULL,
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

  sm <- stanmodels[["MGPCM"]]

  if (is.null(sm)) {
    stop("Stan model MGPCM not found !", call. = FALSE)
  }

  mcmc_env <- setup_mcmc_env(common.method$cores)
  on.exit(restore_mcmc_env(mcmc_env), add = TRUE)

  stan_control <- build_stan_control(algorithm, control.method)

  N <- nrow(response)
  I <- ncol(response)

  response <- istem_prepare_response(response, binary = FALSE)
  length.poly <- istem_prepare_length_poly(response, length.poly)
  max_poly <- max(length.poly)

  Q.matrix <- istem_prepare_01_q(I, D, Q.matrix, triangular = TRUE)

  corr_df <- free_corr_npar(D)

  # Convert response to long format
  data.long <- data.frame(
    y   = as.vector(response),
    pid = rep(1:N, times = I),
    iid = rep(1:I, each = N)
  )

  theta_mu_raw <- get_ctrl("theta.mu", rep(0, D), control)
  theta_mu     <- array(theta_mu_raw, dim = D)

  stan.data <- list(
    N_obs       = nrow(data.long),
    N           = N,
    I           = I,
    D           = D,
    max_poly    = max_poly,
    y           = data.long$y,
    pid         = data.long$pid,
    iid         = data.long$iid,
    length_poly = as.integer(length.poly),
    Q_matrix    = Q.matrix,
    a_mu        = get_ctrl("a.mu", 0.25, control),
    a_sigma     = get_ctrl("a.sigma", 0.25, control),
    d_mu        = get_ctrl("d.mu", 0.0, control),
    d_sigma     = get_ctrl("d.sigma", 1.0, control),
    theta_mu    = theta_mu
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
      show_messages = common.method$vis
    )
  )

  MCMC.obj <- rstan::extract(stan.obj, permuted = TRUE)
  stan.sum <- rstan::summary(stan.obj)$summary

  # ---- Extract theta ----
  theta_est <- extract_theta_stan(MCMC.obj, stan.sum, N, D, rownames(response))

  # ---- Build par matrix: I x (D + max_poly) ----
  # Columns: a1 ... aD, d0, d1, ..., d_{max_poly-1}
  par      <- matrix(NA_real_, nrow = I, ncol = D + max_poly)
  par.se   <- matrix(NA_real_, nrow = I, ncol = D + max_poly)
  par.Rhat <- matrix(NA_real_, nrow = I, ncol = D + max_poly)

  # Extract a parameters
  a_arr  <- ensure_3d(MCMC.obj$a)
  a      <- apply(a_arr, c(2, 3), mean)
  a.se   <- apply(a_arr, c(2, 3), sd)

  par[, 1:D]      <- a
  par.se[, 1:D]   <- a.se
  par.Rhat[, 1:D] <- extract_rhat_matrix(stan.sum, "a", I, D)

  # d0 = 0 (fixed for all items)
  par[, D + 1]      <- 0
  par.se[, D + 1]   <- NA_real_
  par.Rhat[, D + 1] <- NA_real_

  # Extract free d parameters and map back to par matrix
  d_free      <- apply(MCMC.obj$d_free, 2, mean)
  d_free.se   <- apply(MCMC.obj$d_free, 2, sd)
  d_free.idx  <- grep("^d_free\\[", rownames(stan.sum))

  d_pos <- 1  # position in d_free vector
  for (i in 1:I) {
    Ki <- length.poly[i]
    if (Ki > 1) {
      n_d_i <- Ki - 1
      cols <- (D + 2):(D + Ki)
      rows_idx <- d_pos:(d_pos + n_d_i - 1)
      par[i, cols]      <- d_free[rows_idx]
      par.se[i, cols]   <- d_free.se[rows_idx]
      par.Rhat[i, cols] <- stan.sum[d_free.idx[rows_idx], "Rhat"]
      d_pos <- d_pos + n_d_i
    }
  }

  par.free <- mgpcm_par_free_mask(Q.matrix, length.poly, max_poly)
  colnames(par) <- colnames(par.se) <- colnames(par.Rhat) <- colnames(par.free) <-
    c(paste0("a", 1:D), paste0("d", 0:(max_poly - 1)))
  rownames(par) <- rownames(par.se) <- rownames(par.Rhat) <- rownames(par.free) <-
    colnames(response)
  par.se[!par.free] <- NA_real_
  par.Rhat[!par.free] <- NA_real_
  npar <- sum(par.free) + corr_df

  # ---- Extract Corr ----
  Corr_est <- extract_corr_stan(MCMC.obj, stan.sum, D)

  # ---- Assemble results ----
  par_est   <- list(est = par, se = par.se, Rhat = par.Rhat, free = par.free)

  results <- list(
    npar     = npar,
    method   = "stan",
    theta    = theta_est,
    par      = par_est,
    Corr     = Corr_est,
    stan.obj = stan.obj,
    MCMC.obj = MCMC.obj,
    log_lik = MCMC.obj$log_lik,
    Q.matrix    = Q.matrix,
    length.poly = length.poly,
    call        = call,
    arguments = list(
      response  = response,
      D         = D,
      Q.matrix  = Q.matrix,
      length.poly = length.poly,
      method    = "stan",
      chains    = chains,
      iter      = iter,
      warmup    = warmup,
      thin      = thin,
      init      = init,
      algorithm = algorithm,
      cores     = common.method$cores,
      vis       = common.method$vis,
      seed      = common.method$seed,
      control.model = control,
      control.method = fit_effective_method_control(
        control.method, method = stan.method, common = common.method
      )
    )
  )

  class(results) <- "MGPCM"

  L <- get_ctrl("L", NULL, control)
  results$logLik <- logLik.MGPCM(results, L=L)

  return(results)
}
