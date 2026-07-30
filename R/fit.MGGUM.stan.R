#' Fit MGGUM Model via Stan
#'
#' @param response An N x I matrix of integer responses (0, 1, ..., K_i - 1)
#' @param D Number of latent dimensions
#' @param Q.matrix An I x D matrix with values -1, 0, or 1.
#'   \itemize{
#'     \item \code{1}: active dimension with \code{delta} constrained/regularized on the positive side
#'     \item \code{-1}: active dimension with \code{delta} constrained/regularized on the negative side
#'     \item \code{0}: \code{a} and \code{delta} fixed to 0 (inactive dimension)
#'   }
#'   The sign controls the side of the item-location prior/support, not the
#'   sign of discrimination; \code{a} is nonnegative whenever \code{|Q| = 1}.
#'   If NULL, builds a matrix of all 1's (default: positive-side locations).
#' @param length.poly Optional scalar or length-I integer vector giving the
#'   intended number of categories per item. If \code{NULL}, inferred as
#'   \code{max(response[, i]) + 1} for each item.
#' @param control.model Optional named list of model-level controls and
#'   hyperparameters.
#' @param control.method Optional named list of method-specific controls. For
#'   Stan this includes \code{cores}, \code{vis}, \code{seed}, \code{chains},
#'   \code{iter}, \code{warmup}, \code{thin}, \code{init}, \code{algorithm},
#'   and Stan sampler controls.
#'
#' @noRd
fit.MGGUM.stan <- function(response, D = 2, Q.matrix = NULL, length.poly = NULL,
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

  sm <- stanmodels[["MGGUM"]]

  if (is.null(sm)) {
    stop("Stan model MGGUM not found !", call. = FALSE)
  }

  mcmc_env <- setup_mcmc_env(common.method$cores)
  on.exit(restore_mcmc_env(mcmc_env), add = TRUE)

  stan_control <- build_stan_control(algorithm, control.method)

  N <- nrow(response)
  I <- ncol(response)

  response <- istem_prepare_response(response, binary = FALSE)
  length.poly <- istem_prepare_length_poly(response, length.poly)
  max_poly <- max(length.poly)

  Q.matrix <- istem_prepare_signed_q(I, D, Q.matrix)

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
    N_obs           = nrow(data.long),
    N               = N,
    I               = I,
    D               = D,
    max_poly        = max_poly,
    y               = data.long$y,
    pid             = data.long$pid,
    iid             = data.long$iid,
    length_poly     = as.integer(length.poly),
    Q_matrix        = Q.matrix,
    a_mu            = get_ctrl("a.mu", 0, control),
    a_sigma         = get_ctrl("a.sigma", 0.5, control),
    delta_pos_mu    = get_ctrl("delta.pos.mu", 1, control),
    delta_pos_sigma = get_ctrl("delta.pos.sigma", 0.5, control),
    delta_neg_mu    = get_ctrl("delta.neg.mu", -1, control),
    delta_neg_sigma = get_ctrl("delta.neg.sigma", 0.5, control),
    tau_mu          = get_ctrl("tau.mu", 0, control),
    tau_sigma       = get_ctrl("tau.sigma", 0.5, control),
    theta_mu        = theta_mu
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

  # ---- Build par matrix: I x (D + D + max_poly) ----
  # Columns: a1..aD, delta1..deltaD, tau0(=0), tau1, ..., tau_{max_poly-1}
  par      <- matrix(NA_real_, nrow = I, ncol = D + D + max_poly)
  par.se   <- matrix(NA_real_, nrow = I, ncol = D + D + max_poly)
  par.Rhat <- matrix(NA_real_, nrow = I, ncol = D + D + max_poly)

  # Extract a parameters
  a_arr <- ensure_3d(MCMC.obj$a)
  a     <- apply(a_arr, c(2, 3), mean)
  a.se  <- apply(a_arr, c(2, 3), sd)

  par[, 1:D]      <- a
  par.se[, 1:D]   <- a.se
  par.Rhat[, 1:D] <- extract_rhat_matrix(stan.sum, "a", I, D)

  # Extract delta parameters
  delta_arr <- ensure_3d(MCMC.obj$delta)
  delta     <- apply(delta_arr, c(2, 3), mean)
  delta.se  <- apply(delta_arr, c(2, 3), sd)

  par[, (D + 1):(D + D)]      <- delta
  par.se[, (D + 1):(D + D)]   <- delta.se
  par.Rhat[, (D + 1):(D + D)] <- extract_rhat_matrix(stan.sum, "delta", I, D)

  # tau0 = 0 (fixed for all items)
  par[, D + D + 1]      <- 0
  par.se[, D + D + 1]   <- NA_real_
  par.Rhat[, D + D + 1] <- NA_real_

  # Extract tau parameters from MCMC output
  tau_arr <- ensure_3d(MCMC.obj$tau)
  tau     <- apply(tau_arr, c(2, 3), mean)
  tau.se  <- apply(tau_arr, c(2, 3), sd)
  for (i in 1:I) {
    Ki <- length.poly[i]
    if (Ki > 1) {
      tau_cols <- (D + D + 2):(D + D + Ki)
      par[i, tau_cols]      <- tau[i, 2:Ki]
      par.se[i, tau_cols]   <- tau.se[i, 2:Ki]
      # Map back to stan.sum indices: tau[i,2] .. tau[i,Ki]
      tau_idx_i <- grep(paste0("^tau\\[", i, ","), rownames(stan.sum))
      tau_rhat_i <- stan.sum[tau_idx_i, "Rhat"]
      par.Rhat[i, tau_cols] <- tau_rhat_i[2:Ki]
    }
  }

  par.free <- mggum_par_free_mask(Q.matrix, length.poly, max_poly)
  colnames(par) <- colnames(par.se) <- colnames(par.Rhat) <- colnames(par.free) <-
    c(paste0("a", 1:D),
      paste0("delta", 1:D),
      paste0("tau", 0:(max_poly - 1)))
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
      method    = "stan",
      length.poly = length.poly,
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

  class(results) <- "MGGUM"

  L <- get_ctrl("L", NULL, control)
  results$logLik <- logLik.MGGUM(results, L=L)

  return(results)
}
