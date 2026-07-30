#' Fit MIRT Model via Stan
#'
#' @param response An N x I matrix of binary responses (0/1)
#' @param model Model type: \code{"m1pl"}, \code{"m2pl"}, \code{"m3pl"}, or \code{"m4pl"}
#' @param D Number of latent dimensions
#' @param Q.matrix An I x D 0/1 matrix indicating which discrimination parameters
#'   are freely estimated (1) and which are fixed to 0 (default: NULL, builds
#'   a triangular identification structure).
#' @param chains Number of MCMC chains
#' @param iter Number of iterations per chain
#' @param warmup Number of warmup iterations per chain
#' @param thin Thinning interval
#' @param init Initial values (0 = random)
#' @param algorithm Sampling algorithm: \code{"HMC"} (default), \code{"NUTS"},
#'   or \code{"Fixed_param"}.
#' @param control.model A named list of model-level controls and
#'   hyperparameters:
#'   \itemize{
#'     \item \code{a.mu}, \code{a.sigma}: lognormal prior for a parameters
#'     \item \code{b.mu}, \code{b.sigma}: normal prior for item difficulties
#'     \item \code{c.mu}, \code{c.sigma}: uniform support for c parameters
#'     \item \code{d.mu}, \code{d.sigma}: uniform support for d parameters
#'     \item \code{theta.mu}: prior mean vector for theta
#'     \item \code{L}: quadrature grid size used by \code{logLik.MIRT()}
#'   }
#' @param control.method A named list of Stan-specific controls:
#'   \itemize{
#'     \item \code{cores}, \code{vis}, \code{seed}: execution controls
#'     \item \code{chains}, \code{iter}, \code{warmup}, \code{thin},
#'           \code{init}, \code{algorithm}: Stan execution settings
#'     \item \code{adapt_delta}, \code{max_treedepth}, \code{stepsize},
#'           \code{int_time}, \code{metric}: Stan algorithm control parameters
#'   }
#'
#' @noRd
fit.MIRT.stan <- function(response, model = "m2pl", D = 2, Q.matrix = NULL,
                          control.model = NULL, control.method = NULL,
                          .call = NULL) {

  call <- if (is.null(.call)) match.call() else .call

  control <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  stan.method <- fc_stan_method_control(control.method)
  algorithm <- stan.method$algorithm

  model <- tolower(model)

  model_type <- switch(model,
    "m1pl" = 1L,
    "m2pl" = 2L,
    "m3pl" = 3L,
    "m4pl" = 4L,
    stop("Unknown MIRT model type: ", model, call. = FALSE)
  )

  sm <- stanmodels[["MIRT"]]

  if (is.null(sm)) {
    stop("Stan model ", model, " not found !", call. = FALSE)
  }

  mcmc_env <- setup_mcmc_env(common.method$cores)
  on.exit(restore_mcmc_env(mcmc_env), add = TRUE)

  stan_control <- build_stan_control(algorithm, control.method)

  N <- nrow(response)
  I <- ncol(response)

  response <- istem_prepare_response(response, binary = TRUE)
  Q.matrix <- istem_prepare_01_q(I, D, Q.matrix, triangular = TRUE)

  corr_df <- free_corr_npar(D)

  data.long <- data.frame(
    y   = as.vector(response),
    pid = rep(1:N, times = I),
    iid = rep(1:I, each = N)
  )

  theta_mu_raw <- get_ctrl("theta.mu", rep(0, D), control)
  theta_mu     <- array(theta_mu_raw, dim = D)

  stan.data <- list(
    N_obs      = nrow(data.long),
    N          = N,
    I          = I,
    D          = D,
    model_type = model_type,
    Q_matrix   = Q.matrix,
    y          = data.long$y,
    pid        = data.long$pid,
    iid        = data.long$iid,
    a_mu       = get_ctrl("a.mu", 0.25, control),
    a_sigma    = get_ctrl("a.sigma", 0.25, control),
    b_mu       = get_ctrl("b.mu", 0.0, control),
    b_sigma    = get_ctrl("b.sigma", 1.0, control),
    c_mu       = get_ctrl("c.mu", 0.0, control),
    c_sigma    = get_ctrl("c.sigma", 0.35, control),
    d_mu       = get_ctrl("d.mu", 0.65, control),
    d_sigma    = get_ctrl("d.sigma", 1.0, control),
    theta_mu   = theta_mu
  )

  set.seed(common.method$seed)

  stan.obj <- suppressWarnings(
    rstan::sampling(
      object        = sm,
      data          = stan.data,
      chains        = stan.method$chains,
      iter          = stan.method$iter,
      warmup        = stan.method$warmup,
      thin          = stan.method$thin,
      init          = stan.method$init,
      algorithm     = algorithm,
      cores         = mcmc_env$cores,
      control       = stan_control,
      verbose       = common.method$vis,
      show_messages = common.method$vis
    )
  )

  MCMC.obj <- rstan::extract(stan.obj, permuted = TRUE)
  stan.sum <- rstan::summary(stan.obj)$summary

  theta_est <- extract_theta_stan(MCMC.obj, stan.sum, N, D, rownames(response))

  b    <- apply(MCMC.obj$b, 2, mean)
  b.se <- apply(MCMC.obj$b, 2, sd)
  b.idx <- grep("^b\\[", rownames(stan.sum))

  par    <- matrix(NA_real_, nrow = I, ncol = D + 3)
  par.se <- matrix(NA_real_, nrow = I, ncol = D + 3)
  par.Rhat <- matrix(NA_real_, nrow = I, ncol = D + 3)

  if (model %in% c("m2pl", "m3pl", "m4pl")) {
    a_arr  <- ensure_3d(MCMC.obj$a)
    a      <- apply(a_arr, c(2, 3), mean)
    a.se   <- apply(a_arr, c(2, 3), sd)

    par[, 1:D]        <- a
    par.se[, 1:D]     <- a.se
    par.Rhat[, 1:D]   <- extract_rhat_matrix(stan.sum, "a", I, D)

  } else if (model == "m1pl") {
    par[, 1:D]      <- 1
    par.se[, 1:D]   <- NA_real_
    par.Rhat[, 1:D] <- NA_real_
  }

  par[, D + 1]        <- b
  par.se[, D + 1]     <- b.se
  par.Rhat[, D + 1]   <- stan.sum[b.idx, "Rhat"]

  if (model %in% c("m3pl", "m4pl")) {
    c_est    <- apply(MCMC.obj$c, 2, mean)
    c_est.se <- apply(MCMC.obj$c, 2, sd)
    c.idx <- grep("^c\\[", rownames(stan.sum))
    par[, D + 2]      <- c_est
    par.se[, D + 2]   <- c_est.se
    par.Rhat[, D + 2] <- stan.sum[c.idx, "Rhat"]
  } else {
    par[, D + 2]      <- 0
    par.se[, D + 2]   <- NA_real_
    par.Rhat[, D + 2] <- NA_real_
  }

  if (model == "m4pl") {
    d_est    <- apply(MCMC.obj$d, 2, mean)
    d_est.se <- apply(MCMC.obj$d, 2, sd)
    d.idx <- grep("^d\\[", rownames(stan.sum))
    par[, D + 3]      <- d_est
    par.se[, D + 3]   <- d_est.se
    par.Rhat[, D + 3] <- stan.sum[d.idx, "Rhat"]
  } else {
    par[, D + 3]      <- 1
    par.se[, D + 3]   <- NA_real_
    par.Rhat[, D + 3] <- NA_real_
  }

  par.free <- mirt_par_free_mask(model_type, Q.matrix)
  colnames(par) <- colnames(par.se) <- colnames(par.Rhat) <- colnames(par.free) <-
    c(paste0("a", 1:D), "b", "c", "d")
  rownames(par) <- rownames(par.se) <- rownames(par.Rhat) <- rownames(par.free) <-
    colnames(response)
  par.se[!par.free] <- NA_real_
  par.Rhat[!par.free] <- NA_real_
  npar <- sum(par.free) + corr_df

  Corr_est <- extract_corr_stan(MCMC.obj, stan.sum, D)
  par_est  <- list(est = par, se = par.se, Rhat = par.Rhat, free = par.free)

  results <- list(
    npar     = npar,
    method   = "stan",
    theta    = theta_est,
    par      = par_est,
    Corr     = Corr_est,
    stan.obj = stan.obj,
    MCMC.obj = MCMC.obj,
    log_lik = MCMC.obj$log_lik,
    Q.matrix = Q.matrix,
    call = call,
    arguments = list(
      response  = response,
      model     = model,
      D         = D,
      Q.matrix  = Q.matrix,
      method    = "stan",
      cores     = common.method$cores,
      vis       = common.method$vis,
      seed      = common.method$seed,
      control.model = control,
      control.method = fit_effective_method_control(
        control.method, method = stan.method, common = common.method
      )
    )
  )

  class(results) <- "MIRT"

  L <- get_ctrl("L", NULL, control)
  results$logLik <- logLik.MIRT(results, L=L)

  return(results)
}
