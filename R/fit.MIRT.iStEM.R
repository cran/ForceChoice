#' Fit MIRT Model via iStEM
#'
#' iStEM uses the same default item-prior family as the Stan backend in the
#' item MAP update: \code{a ~ lognormal(a.mu, a.sigma)},
#' \code{b ~ normal(b.mu, b.sigma)}, \code{c ~ uniform(c.mu, c.sigma)}, and
#' \code{d ~ uniform(d.mu, d.sigma)}. Set \code{use.prior = FALSE} in
#' \code{control.model} to use the previous unpenalized item update.
#'
#' @noRd
fit.MIRT.iStEM <- function(response, model = "m2pl", D = 2, Q.matrix = NULL,
                           control.model = NULL, control.method = NULL,
                           .call = NULL) {

  call <- if (is.null(.call)) match.call() else .call
  control.model <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  vis <- common.method$vis
  set.seed(common.method$seed)

  response <- as.matrix(response)
  if (anyNA(response) || !all(response %in% c(0, 1))) {
    stop("'response' must be an N x I matrix of binary responses (0/1).", call. = FALSE)
  }
  storage.mode(response) <- "integer"

  model <- tolower(model)
  model_type <- model_type_to_int(model)

  N <- nrow(response)
  I <- ncol(response)
  response.group <- istem_response_groups(response)
  response.fit <- response.group$response
  count.fit <- response.group$count
  D <- as.integer(D)
  if (D < 1L) {
    stop("'D' must be a positive integer.", call. = FALSE)
  }

  q_info <- mirt_istem_prepare_q(I, D, Q.matrix, model)
  Q.matrix <- q_info$Q.matrix

  method <- istem_method_control(control.method)
  method <- istem_apply_grid_control(method, control.model, D)
  istem_check_method_control(method)
  control.model <- istem_effective_grid_model_control(control.model, method)
  M            <- method$M
  B            <- method$B
  burnin.maxitr <- method$burnin.maxitr
  maxitr       <- method$maxitr
  eps1         <- method$eps1
  eps2         <- method$eps2
  frac1        <- method$frac1
  frac2        <- method$frac2
  L            <- method$L
  theta_lower  <- method$theta.lower
  theta_upper  <- method$theta.upper
  optim_maxit  <- method$optim.maxit
  corr_optim_maxit <- method$corr.optim.maxit
  fix.corr     <- method$fix.corr
  estimate.se  <- method$estimate.se
  include_corr <- D > 1L && !fix.corr
  corr_df      <- free_corr_npar(D, include_corr)

  theta_mu_raw <- get_ctrl("theta.mu", rep(0, D), control.model)
  theta_mu <- array(theta_mu_raw, dim = D)
  prior <- mirt_istem_prior_control(control.model)

  bounds <- list(
    a.lower = get_ctrl("a.lower", 1e-4, control.method),
    a.upper = get_ctrl("a.upper", 6, control.method),
    b.lower = get_ctrl("b.lower", -6, control.method),
    b.upper = get_ctrl("b.upper", 6, control.method),
    c.lower = get_ctrl("c.lower", prior$c.lower, control.method),
    c.upper = get_ctrl("c.upper", prior$c.upper, control.method),
    d.lower = get_ctrl("d.lower", prior$d.lower, control.method),
    d.upper = get_ctrl("d.upper", prior$d.upper, control.method)
  )
  mirt_istem_check_bounds(bounds, model_type)

  istem.method <- c(
    list(
      M = M,
      B = B,
      burnin.maxitr = burnin.maxitr,
      maxitr = maxitr,
      eps1 = eps1,
      eps2 = eps2,
      frac1 = frac1,
      frac2 = frac2,
      L = L,
      theta.lower = theta_lower,
      theta.upper = theta_upper,
      optim.maxit = optim_maxit,
      corr.optim.maxit = corr_optim_maxit,
      fix.corr = fix.corr,
      estimate.se = estimate.se
    ),
    bounds
  )

  par <- mirt_istem_initial_par(response.fit, model_type, Q.matrix, bounds, prior)
  theta <- istem_random_theta(response.group$G, D, theta_mu,
                              theta_lower, theta_upper)
  Corr <- istem_random_corr(D, include_corr)

  state <- list(
    response = response.fit,
    response.count = count.fit,
    par = par,
    theta = theta,
    Corr = Corr,
    theta_mu = theta_mu,
    Q.matrix = Q.matrix,
    model_type = model_type,
    optim_maxit = optim_maxit,
    bounds = bounds,
    prior = prior
  )

  run <- istem_run(
    state = state,
    method = istem.method,
    N = N,
    vis = vis,
    label = "MIRT",
    update_parameters = function(state) {
      for (j in seq_len(ncol(state$response))) {
        state$par[j, ] <- mirt_istem_fit_item(
          y = state$response[, j],
          theta = state$theta,
          par_i = state$par[j, ],
          model_type = state$model_type,
          q_i = state$Q.matrix[j, ],
          weight = istem_state_weight(state),
          optim_maxit = state$optim_maxit,
          bounds = state$bounds,
          prior = state$prior
        )
      }
      state
    },
    param_vec = function(state) {
      mirt_istem_param_vec(state$par, state$Corr, state$model_type,
                           state$Q.matrix, include_corr)
    },
    sample_theta = function(state, theta_grid_length, theta_lower, theta_upper) {
      mirt_istem_sample_theta_batch(
        theta = state$theta,
        response = state$response,
        par = state$par,
        Corr = state$Corr,
        theta_mu = state$theta_mu,
        grid_length = theta_grid_length,
        lower = theta_lower,
        upper = theta_upper
      )
    },
    logLik_fun = function(state) {
      mirt_istem_loglik_trace(state, L, theta_lower, theta_upper)
    }
  )

  final <- mirt_istem_param_unpack(
    pv = run$chain.mean,
    pv_se = run$chain.sd,
    pv_rhat = NULL,
    I = I,
    D = D,
    model_type = model_type,
    Q.matrix = Q.matrix,
    Corr = run$state$Corr,
    include_corr = include_corr
  )

  theta.est <- istem_expand_group_matrix(
    run$theta.est, response.group$group, rownames(response))
  theta.se <- istem_expand_group_matrix(
    run$theta.se, response.group$group, rownames(response))
  theta <- istem_named_theta(theta.est, theta.se, response, D)

  par <- final$par
  par.se <- final$par.se
  par.Rhat <- final$par.rhat
  par.free <- mirt_par_free_mask(model_type, Q.matrix)
  colnames(par) <- colnames(par.se) <- colnames(par.Rhat) <- colnames(par.free) <-
    c(paste0("a", seq_len(D)), "b", "c", "d")
  rownames(par) <- rownames(par.se) <- rownames(par.Rhat) <- rownames(par.free) <-
    colnames(response)
  par.se[!par.free] <- NA_real_
  par.Rhat[!par.free] <- NA_real_
  npar <- sum(par.free) + corr_df

  Corr.Rhat <- final$Corr.rhat
  dimnames(Corr.Rhat) <- list(paste0("Dim.", seq_len(D)), paste0("Dim.", seq_len(D)))
  Corr <- list(est = final$Corr, se = final$Corr.se, Rhat = Corr.Rhat)

  results <- list(
    npar = npar,
    method = "iStEM",
    theta = theta,
    par = list(est = par, se = par.se, Rhat = par.Rhat, free = par.free),
    Corr = Corr,
    stan.obj = NULL,
    MCMC.obj = NULL,
    Q.matrix = Q.matrix,
    call = call,
    arguments = list(
      response = response,
      model = model,
      D = D,
      Q.matrix = Q.matrix,
      method = "iStEM",
      cores = common.method$cores,
      vis = vis,
      seed = common.method$seed,
      control.model = control.model,
      control.method = fit_effective_method_control(
        control.method, method = istem.method, common = common.method
      )
    ),
    iStEM = utils::modifyList(run$iStEM, list(prior = prior))
  )

  class(results) <- "MIRT"
  results$logLik <- logLik.MIRT(results, L = method$L,
                                theta.low = method$theta.lower,
                                theta.up = method$theta.upper)

  results
}

mirt_istem_prepare_q <- function(I, D, Q.matrix, model) {
  if (!is.null(Q.matrix)) {
    Q.matrix <- as.matrix(Q.matrix)
    if (ncol(Q.matrix) != D || nrow(Q.matrix) != I) {
      Q.matrix <- NULL
    }
  }

  if (is.null(Q.matrix)) {
    Q.matrix <- matrix(1, I, D)
    if (D > 1L) {
      for (k in seq_len(D)) {
        r <- I - D + k
        for (j in seq_len(D)) {
          if (k > D - j + 1L) {
            Q.matrix[r, j] <- 0
          }
        }
      }
    }
  }

  list(Q.matrix = Q.matrix)
}

mirt_istem_prior_control <- function(control.model) {
  use_prior <- get_ctrl("use.prior", TRUE, control.model)
  if (!is.logical(use_prior) || length(use_prior) != 1L || is.na(use_prior)) {
    stop("'use.prior' in 'control.model' must be TRUE or FALSE.",
         call. = FALSE)
  }

  a_mu <- istem_scalar(get_ctrl("a.mu", 0.25, control.model), "a.mu")
  a_sigma <- istem_scalar(
    get_ctrl("a.sigma", 0.25, control.model),
    "a.sigma",
    positive = TRUE
  )
  b_mu <- istem_scalar(get_ctrl("b.mu", 0.0, control.model), "b.mu")
  b_sigma <- istem_scalar(
    get_ctrl("b.sigma", 1.0, control.model),
    "b.sigma",
    positive = TRUE
  )
  c_lower <- istem_scalar(get_ctrl("c.mu", 0.0, control.model), "c.mu")
  c_upper <- istem_scalar(get_ctrl("c.sigma", 0.35, control.model), "c.sigma")
  d_lower <- istem_scalar(get_ctrl("d.mu", 0.65, control.model), "d.mu")
  d_upper <- istem_scalar(get_ctrl("d.sigma", 1.0, control.model), "d.sigma")

  if (c_lower >= c_upper) {
    stop("'c.mu' must be smaller than 'c.sigma' for the c uniform support.",
         call. = FALSE)
  }
  if (d_lower >= d_upper) {
    stop("'d.mu' must be smaller than 'd.sigma' for the d uniform support.",
         call. = FALSE)
  }

  list(
    use.prior = use_prior,
    a.mu = a_mu,
    a.sigma = a_sigma,
    b.mu = b_mu,
    b.sigma = b_sigma,
    c.mu = c_lower,
    c.sigma = c_upper,
    d.mu = d_lower,
    d.sigma = d_upper,
    c.lower = c_lower,
    c.upper = c_upper,
    d.lower = d_lower,
    d.upper = d_upper
  )
}

mirt_istem_check_bounds <- function(bounds, model_type) {
  for (nm in names(bounds)) {
    if (!is.numeric(bounds[[nm]]) || length(bounds[[nm]]) != 1L ||
        !is.finite(bounds[[nm]])) {
      stop("'", nm, "' in 'control.method' must be a finite numeric scalar.",
           call. = FALSE)
    }
    bounds[[nm]] <- as.numeric(bounds[[nm]])
  }

  if (model_type >= 2L && bounds$a.lower <= 0) {
    stop("'a.lower' in 'control.method' must be positive for iStEM.",
         call. = FALSE)
  }
  if (bounds$a.lower >= bounds$a.upper) {
    stop("'a.lower' must be smaller than 'a.upper'.", call. = FALSE)
  }
  if (bounds$b.lower >= bounds$b.upper) {
    stop("'b.lower' must be smaller than 'b.upper'.", call. = FALSE)
  }
  if (model_type >= 3L && bounds$c.lower >= bounds$c.upper) {
    stop("'c.lower' must be smaller than 'c.upper'.", call. = FALSE)
  }
  if (model_type == 4L && bounds$d.lower >= bounds$d.upper) {
    stop("'d.lower' must be smaller than 'd.upper'.", call. = FALSE)
  }
  invisible(NULL)
}

mirt_istem_initial_par <- function(response, model_type, Q.matrix, bounds,
                                   prior) {
  I <- ncol(response)
  D <- ncol(Q.matrix)
  b <- istem_rnorm_bounded(I, prior$b.mu, prior$b.sigma,
                           bounds$b.lower, bounds$b.upper)

  if (model_type == 1L) {
    a <- matrix(1, I, D)
  } else {
    a <- matrix(0, I, D)
    a[Q.matrix == 1] <- istem_rlnorm_bounded(
      sum(Q.matrix == 1), prior$a.mu, prior$a.sigma,
      bounds$a.lower, bounds$a.upper
    )
  }

  c_par <- if (model_type >= 3L) {
    stats::runif(I, max(bounds$c.lower, prior$c.lower),
                 min(bounds$c.upper, prior$c.upper))
  } else {
    rep(0, I)
  }
  d_par <- if (model_type == 4L) {
    stats::runif(I, max(bounds$d.lower, prior$d.lower),
                 min(bounds$d.upper, prior$d.upper))
  } else {
    rep(1, I)
  }

  cbind(a, b, c_par, d_par)
}

mirt_istem_sample_theta_batch <- function(theta, response, par, Corr, theta_mu,
                                          grid_length, lower, upper) {
  chol <- istem_corr_chol(Corr)
  cpp_gibbs_mirt_theta(
    theta_ = theta,
    par_ = par,
    response_ = response,
    theta_mu_ = theta_mu,
    chol_corr_ = chol$chol,
    log_diag_sum = chol$log_diag_sum,
    step = grid_length,
    lower = lower,
    upper = upper
  )
}

mirt_istem_fit_item <- function(y, theta, par_i, model_type, q_i, weight,
                                optim_maxit, bounds, prior) {
  D <- ncol(theta)
  active <- q_i == 1

  if (model_type == 1L) {
    x0 <- par_i[D + 1L]
    lower <- bounds$b.lower
    upper <- bounds$b.upper
  } else {
    x0 <- c(par_i[seq_len(D)][active], par_i[D + 1L])
    lower <- c(rep(bounds$a.lower, sum(active)), bounds$b.lower)
    upper <- c(rep(bounds$a.upper, sum(active)), bounds$b.upper)

    if (model_type >= 3L) {
      x0 <- c(x0, par_i[D + 2L])
      lower <- c(lower, bounds$c.lower)
      upper <- c(upper, bounds$c.upper)
    }
    if (model_type == 4L) {
      x0 <- c(x0, par_i[D + 3L])
      lower <- c(lower, bounds$d.lower)
      upper <- c(upper, bounds$d.upper)
    }
  }

  obj <- function(x) {
    item <- mirt_istem_x_to_item(x, par_i, model_type, active, D)
    ll <- cpp_mirt_item_loglik_weighted(
      theta = theta,
      a = item$a,
      b = item$b,
      c = item$c,
      upper = item$d,
      response = as.integer(y),
      weight = weight
    )
    log_prior <- mirt_istem_item_log_prior(
      item = item,
      model_type = model_type,
      active = active,
      prior = prior
    )
    if (!is.finite(log_prior)) {
      return(Inf)
    }
    -ll - log_prior
  }

  opt <- tryCatch(
    stats::optim(
      par = x0,
      fn = obj,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = optim_maxit)
    ),
    error = function(e) NULL
  )

  if (is.null(opt) || !is.finite(opt$value)) {
    return(par_i)
  }

  item <- mirt_istem_x_to_item(opt$par, par_i, model_type, active, D)
  c(item$a, item$b, item$c, item$d)
}

mirt_istem_item_log_prior <- function(item, model_type, active, prior) {
  if (!isTRUE(prior$use.prior)) {
    return(0)
  }

  lp <- stats::dnorm(item$b, mean = prior$b.mu, sd = prior$b.sigma, log = TRUE)
  if (model_type >= 2L) {
    a_active <- item$a[active]
    if (any(a_active <= 0)) {
      return(-Inf)
    }
    lp <- lp + sum(stats::dlnorm(
      a_active,
      meanlog = prior$a.mu,
      sdlog = prior$a.sigma,
      log = TRUE
    ))
  }
  if (model_type >= 3L &&
      (item$c < prior$c.lower || item$c > prior$c.upper)) {
    return(-Inf)
  }
  if (model_type == 4L &&
      (item$d < prior$d.lower || item$d > prior$d.upper)) {
    return(-Inf)
  }

  lp
}

mirt_istem_x_to_item <- function(x, par_i, model_type, active, D) {
  idx <- 1L
  if (model_type == 1L) {
    a <- rep(1, D)
    b <- x[1L]
    c_par <- 0
    d_par <- 1
  } else {
    a <- rep(0, D)
    n_active <- sum(active)
    if (n_active > 0L) {
      a[active] <- x[idx:(idx + n_active - 1L)]
      idx <- idx + n_active
    }
    b <- x[idx]
    idx <- idx + 1L
    c_par <- if (model_type >= 3L) x[idx] else 0
    if (model_type >= 3L) {
      idx <- idx + 1L
    }
    d_par <- if (model_type == 4L) x[idx] else 1
  }
  list(a = a, b = b, c = c_par, d = d_par)
}

mirt_istem_param_vec <- function(par, Corr, model_type, Q.matrix, include_corr) {
  D <- ncol(Q.matrix)
  pv <- if (model_type == 1L) {
    par[, D + 1L]
  } else {
    c(par[, seq_len(D), drop = FALSE][Q.matrix == 1], par[, D + 1L])
  }

  if (model_type >= 3L) {
    pv <- c(pv, par[, D + 2L])
  }
  if (model_type == 4L) {
    pv <- c(pv, par[, D + 3L])
  }
  if (include_corr) {
    pv <- c(pv, Corr[lower.tri(Corr)])
  }
  pv
}

mirt_istem_param_unpack <- function(pv, pv_se, pv_rhat = NULL, I, D, model_type,
                                    Q.matrix, Corr, include_corr) {
  idx <- 1L
  par <- matrix(NA_real_, I, D + 3L)
  par.se <- matrix(NA_real_, I, D + 3L)
  par.rhat <- matrix(NA_real_, I, D + 3L)
  has_rhat <- !is.null(pv_rhat)

  if (model_type == 1L) {
    par[, seq_len(D)] <- 1
    par.se[, seq_len(D)] <- NA_real_
    par[, D + 1L] <- pv[idx:(idx + I - 1L)]
    par.se[, D + 1L] <- pv_se[idx:(idx + I - 1L)]
    if (has_rhat) par.rhat[, D + 1L] <- pv_rhat[idx:(idx + I - 1L)]
    idx <- idx + I
  } else {
    a <- matrix(0, I, D)
    a.se <- matrix(NA_real_, I, D)
    a.rhat <- matrix(NA_real_, I, D)
    n_active <- sum(Q.matrix == 1)
    a[Q.matrix == 1] <- pv[idx:(idx + n_active - 1L)]
    a.se[Q.matrix == 1] <- pv_se[idx:(idx + n_active - 1L)]
    if (has_rhat) a.rhat[Q.matrix == 1] <- pv_rhat[idx:(idx + n_active - 1L)]
    idx <- idx + n_active
    par[, seq_len(D)] <- a
    par.se[, seq_len(D)] <- a.se
    par.rhat[, seq_len(D)] <- a.rhat
    par[, D + 1L] <- pv[idx:(idx + I - 1L)]
    par.se[, D + 1L] <- pv_se[idx:(idx + I - 1L)]
    if (has_rhat) par.rhat[, D + 1L] <- pv_rhat[idx:(idx + I - 1L)]
    idx <- idx + I
  }

  if (model_type >= 3L) {
    par[, D + 2L] <- pv[idx:(idx + I - 1L)]
    par.se[, D + 2L] <- pv_se[idx:(idx + I - 1L)]
    if (has_rhat) par.rhat[, D + 2L] <- pv_rhat[idx:(idx + I - 1L)]
    idx <- idx + I
  } else {
    par[, D + 2L] <- 0
    par.se[, D + 2L] <- NA_real_
  }

  if (model_type == 4L) {
    par[, D + 3L] <- pv[idx:(idx + I - 1L)]
    par.se[, D + 3L] <- pv_se[idx:(idx + I - 1L)]
    if (has_rhat) par.rhat[, D + 3L] <- pv_rhat[idx:(idx + I - 1L)]
    idx <- idx + I
  } else {
    par[, D + 3L] <- 1
    par.se[, D + 3L] <- NA_real_
  }

  Corr.se <- matrix(NA_real_, D, D)
  if (D == 1L) {
    Corr <- matrix(1, 1, 1)
    Corr.se <- matrix(0, 1, 1)
  } else if (include_corr) {
    n_corr <- D * (D - 1L) / 2L
    corr_vals <- pv[idx:(idx + n_corr - 1L)]
    corr_se <- pv_se[idx:(idx + n_corr - 1L)]
    Corr <- diag(D)
    Corr[lower.tri(Corr)] <- corr_vals
    Corr <- t(Corr)
    Corr[lower.tri(Corr)] <- corr_vals
    Corr.se[lower.tri(Corr.se)] <- corr_se
    Corr.se <- t(Corr.se)
    Corr.se[lower.tri(Corr.se)] <- corr_se
    diag(Corr.se) <- 0
    Corr <- istem_safe_corr(Corr)
  } else {
    Corr.se[,] <- NA_real_
    diag(Corr.se) <- 0
  }
  Corr.rhat <- istem_unpack_corr_rhat(pv_rhat, idx, D, include_corr)

  list(par = par, par.se = par.se, par.rhat = par.rhat,
       Corr = Corr, Corr.se = Corr.se, Corr.rhat = Corr.rhat)
}
