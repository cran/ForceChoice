# ---- MCMC environment & algorithm helpers -------------------------------------

#' Setup single-threaded MCMC environment
#'
#' Ensures each MCMC chain runs strictly on a single OS thread by setting
#' environment variables that control Stan-internal threading and BLAS-level
#' parallelism.  The validated core count is returned inside the result list.
#'
#' @param cores Desired number of parallel chains (one chain per core).
#' @return A named list with old environment values and the validated
#'   \code{cores} value.  Pass this to \code{restore_mcmc_env()} to restore
#'   the previous thread-control state.
#' @keywords internal
setup_mcmc_env <- function(cores) {
  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1L) {
    stop("'cores' must be a positive integer, got ", cores, call. = FALSE)
  }

  max_cores <- parallel::detectCores(logical = TRUE)
  if (cores > max_cores) {
    warning(sprintf(
      "'cores' (%d) exceeds available logical cores (%d). Reducing to %d.",
      cores, max_cores, max_cores
    ))
    cores <- max_cores
  }

  old <- list(
    STAN_NUM_THREADS     = Sys.getenv("STAN_NUM_THREADS",     unset = NA_character_),
    OMP_NUM_THREADS      = Sys.getenv("OMP_NUM_THREADS",      unset = NA_character_),
    MKL_NUM_THREADS      = Sys.getenv("MKL_NUM_THREADS",      unset = NA_character_),
    OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA_character_),
    cores                = cores
  )

  Sys.setenv(STAN_NUM_THREADS     = "1")
  Sys.setenv(OMP_NUM_THREADS      = "1")
  Sys.setenv(MKL_NUM_THREADS      = "1")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")

  return(old)
}


#' Restore MCMC threading environment to its previous state
#'
#' @param old A list returned by \code{setup_mcmc_env()}.
#' @return Invisibly \code{NULL}.
#' @keywords internal
restore_mcmc_env <- function(old) {
  restore_one <- function(name) {
    val <- old[[name]]
    if (is.na(val)) {
      Sys.unsetenv(name)
    } else {
      args <- list(val)
      names(args) <- name
      do.call(Sys.setenv, args)
    }
  }
  restore_one("STAN_NUM_THREADS")
  restore_one("OMP_NUM_THREADS")
  restore_one("MKL_NUM_THREADS")
  restore_one("OPENBLAS_NUM_THREADS")
  invisible(NULL)
}


#' Build Stan sampling-control list from user-provided \code{control.method}
#'
#' This helper extracts Stan sampler tuning knobs (e.g. \code{adapt_delta},
#' \code{max_treedepth}) from method-specific controls so they can be passed
#' safely to \code{rstan::sampling()}.
#'
#' @param algorithm Character scalar: \code{"NUTS"}, \code{"HMC"}, or
#'   \code{"Fixed_param"}.
#' @param control.method A named list of Stan method controls.
#' @return A named list of Stan sampling-control parameters (possibly empty).
#' @keywords internal
build_stan_control <- function(algorithm, control.method) {
  if (is.null(control.method)) control.method <- list()

  stan_ctrl_names <- c(
    "adapt_delta", "stepsize", "max_treedepth",
    "int_time", "metric", "adapt_engaged",
    "adapt_init_buffer", "adapt_term_buffer", "adapt_window"
  )

  stan_ctrl <- list(
    adapt_delta   = 0.95,
    max_treedepth = 12L
  )
  for (nm in stan_ctrl_names) {
    if (!is.null(control.method[[nm]])) {
      stan_ctrl[[nm]] <- control.method[[nm]]
    }
  }

  if (algorithm == "Fixed_param") {
    stan_ctrl$adapt_engaged <- FALSE
  }

  return(stan_ctrl)
}


# ---- Stan MCMC post-processing helpers ---------------------------------------

#' Extract theta posterior summary from Stan MCMC output
#'
#' @param MCMC.obj Output of \code{rstan::extract()}.
#' @param stan.sum Output of \code{rstan::summary()$summary}.
#' @param N Number of persons.
#' @param D Number of dimensions.
#' @param row_names Row names for the output matrices (typically person IDs).
#' @return A list with components \code{est}, \code{se}, and \code{Rhat}.
#' @keywords internal
extract_theta_stan <- function(MCMC.obj, stan.sum, N, D, row_names = NULL) {
  theta_arr  <- ensure_3d(MCMC.obj$theta)
  theta      <- apply(theta_arr, c(2, 3), mean)
  theta.se   <- apply(theta_arr, c(2, 3), sd)
  theta.Rhat <- extract_rhat_matrix(stan.sum, "theta", N, D)
  dn <- list(row_names, paste0("Dim.", seq_len(D)))
  dimnames(theta) <- dimnames(theta.se) <- dimnames(theta.Rhat) <- dn
  list(est = theta, se = theta.se, Rhat = theta.Rhat)
}

#' Extract correlation posterior summary from Stan MCMC output
#'
#' @param MCMC.obj Output of \code{rstan::extract()}.
#' @param stan.sum Output of \code{rstan::summary()$summary}.
#' @param D Number of dimensions.
#' @return A list with components \code{est}, \code{se}, and \code{Rhat}.
#' @keywords internal
extract_corr_stan <- function(MCMC.obj, stan.sum, D) {
  Corr_arr  <- ensure_3d(MCMC.obj$Corr)
  Corr      <- apply(Corr_arr, c(2, 3), mean)
  Corr.se   <- apply(Corr_arr, c(2, 3), sd)
  Corr.Rhat <- extract_rhat_matrix(stan.sum, "Corr", D, D)
  diag(Corr.Rhat) <- NA_real_
  dn <- list(paste0("Dim.", seq_len(D)), paste0("Dim.", seq_len(D)))
  dimnames(Corr) <- dimnames(Corr.se) <- dimnames(Corr.Rhat) <- dn
  list(est = Corr, se = Corr.se, Rhat = Corr.Rhat)
}


# ---- Quadrature grid helpers for logLik ---------------------------------------

#' Multivariate normal quadrature grid with Cholesky-stabilised weights
#'
#' Used by all continuous-trait \code{logLik()} methods.
#'
#' @param Corr D x D correlation matrix (can be \code{NULL} for identity).
#' @param theta.mu Numeric vector of length D (latent mean).
#' @param D Number of dimensions.
#' @param L Number of quadrature points per dimension, or \code{NULL} to
#'   auto-select from a built-in table.
#' @param theta.low,theta.up Quadrature range.
#' @return A list with \code{theta.norm} (grid matrix) and \code{pi}
#'   (normalised prior-weight column vector).
#' @keywords internal
make_quadrature_grid_mvn <- function(Corr, theta.mu, D, L = NULL,
                                     theta.low = -6, theta.up = 6) {
  if (is.null(L)) {
    L <- switch(as.character(D),
                "1" = 61, "2" = 31, "3" = 15, "4" = 9, "5" = 7, 3)
  }

  theta.norm.list <- vector("list", D)
  theta.axis <- seq(theta.low, theta.up, length.out = L)
  for (d in seq_len(D)) {
    theta.norm.list[[d]] <- theta.axis
  }
  theta.norm <- as.matrix(do.call(expand.grid, theta.norm.list))

  if (is.null(Corr) || anyNA(Corr)) Corr <- diag(D)
  Corr <- as.matrix(Corr)
  diag(Corr) <- 1

  logspace_sum <- function(x) {
    m <- max(x)
    if (!is.finite(m)) return(m)
    m + log(sum(exp(x - m)))
  }

  if (D == 1L) {
    log.pi.raw <- stats::dnorm(theta.norm[, 1L], mean = theta.mu[1L],
                               sd = 1, log = TRUE)
  } else {
    chol.ok <- FALSE
    jitter <- 0
    while (!chol.ok && jitter <= 1e-4) {
      chol.try <- try(chol(Corr + diag(jitter, D)), silent = TRUE)
      if (!inherits(chol.try, "try-error")) {
        chol.ok <- TRUE
      } else {
        jitter <- if (jitter == 0) 1e-10 else jitter * 10
      }
    }
    if (!chol.ok) {
      stop("Estimated 'Corr' is not positive definite enough for quadrature.", call. = FALSE)
    }
    centered.theta <- sweep(theta.norm, 2, theta.mu, "-")
    z.theta <- forwardsolve(t(chol.try), t(centered.theta))
    log.det <- 2 * sum(log(diag(chol.try)))
    quad <- colSums(z.theta^2)
    log.pi.raw <- -0.5 * (D * log(2 * pi) + log.det + quad)
  }
  log.pi.raw <- log.pi.raw - logspace_sum(log.pi.raw)
  pi <- matrix(exp(log.pi.raw), ncol = 1)

  list(theta.norm = theta.norm, pi = pi, L = L)
}

#' Apply response-pattern frequencies to grouped log-likelihoods
#'
#' @param lik A grouped likelihood list containing \code{log_marginal}.
#' @param count Frequency of each unique response pattern.
#' @return \code{lik} with aggregate \code{logLik} and
#'   \code{response.count} components.
#' @noRd
loglik_apply_group_count <- function(lik, count) {
  count <- as.numeric(count)
  if (is.null(lik$log_marginal) || length(lik$log_marginal) != length(count)) {
    stop("'lik' must contain one marginal log-likelihood per response group.",
         call. = FALSE)
  }
  lik$logLik <- sum(lik$log_marginal * count)
  lik$response.count <- count
  lik
}

#' Attach standard attributes to a logLik result vector
#'
#' @param results A named numeric vector containing the log-likelihood.
#' @param N Number of observations.
#' @param df Number of free parameters.
#' @param lik A list returned by a \code{cpp_loglik_*()} function.
#' @param theta.norm Quadrature grid matrix.
#' @param prob Model-implied probability matrix.
#' @param pi Quadrature or latent-class weights.
#' @param L Number of quadrature points per dimension.
#' @param theta.low,theta.up Quadrature range.
#' @param ... Additional named attributes.
#' @return \code{results} with class \code{"logLik"} and standard
#'   likelihood attributes.
#' @noRd
set_loglik_attributes <- function(results, N, df, lik, theta.norm,
                                   prob, pi, L, theta.low, theta.up, ...) {
  class(results) <- "logLik"
  attr(results, "nobs")        <- N
  attr(results, "df")          <- df
  attr(results, "L.theta.Xi")  <- lik$L.theta.Xi
  attr(results, "P.theta.Xi")  <- lik$P.theta.Xi
  attr(results, "theta.norm")  <- theta.norm
  attr(results, "prob")        <- prob
  attr(results, "pi")          <- pi
  attr(results, "L")           <- L
  attr(results, "theta.low")   <- theta.low
  attr(results, "theta.up")    <- theta.up

  dots <- list(...)
  if (length(dots) > 0L) {
    for (nm in names(dots)) {
      attr(results, nm) <- dots[[nm]]
    }
  }

  results
}
