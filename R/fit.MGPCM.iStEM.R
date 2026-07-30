################################# MGPCM iStEM #################################

fit.MGPCM.iStEM <- function(response, D = 2, Q.matrix = NULL,
                            length.poly = NULL,
                            control.model = NULL, control.method = NULL,
                            .call = NULL) {

  call <- if (is.null(.call)) match.call() else .call
  control.model <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  vis <- common.method$vis
  method <- istem_method_control(control.method)
  set.seed(common.method$seed)

  response <- as.matrix(response)
  if (anyNA(response) || any(response < 0) ||
      any(response != floor(response))) {
    stop("'response' must be an N x I matrix of integer responses starting at 0.",
         call. = FALSE)
  }
  storage.mode(response) <- "integer"

  N <- nrow(response)
  I <- ncol(response)
  response.group <- istem_response_groups(response)
  response.fit <- response.group$response
  count.fit <- response.group$count
  D <- as.integer(D)
  if (D < 1L) {
    stop("'D' must be a positive integer.", call. = FALSE)
  }
  method <- istem_apply_grid_control(method, control.model, D)
  istem_check_method_control(method)
  control.model <- istem_effective_grid_model_control(control.model, method)

  length.poly <- istem_prepare_length_poly(response, length.poly)
  max_poly <- max(length.poly)
  Q.matrix <- istem_prepare_01_q(I, D, Q.matrix, triangular = TRUE)

  prior <- mgpcm_istem_prior_control(control.model)
  bounds <- list(
    a.lower = get_ctrl("a.lower", 1e-4, control.method),
    a.upper = get_ctrl("a.upper", 6, control.method),
    d.lower = get_ctrl("d.lower", -8, control.method),
    d.upper = get_ctrl("d.upper", 8, control.method)
  )
  mgpcm_istem_check_bounds(bounds)

  theta_mu <- array(get_ctrl("theta.mu", rep(0, D), control.model), dim = D)
  include_corr <- D > 1L && !method$fix.corr
  corr_df <- free_corr_npar(D, include_corr)
  par <- mgpcm_istem_initial_par(response.fit, Q.matrix, length.poly, bounds,
                                 prior)

  state <- list(
    response = response.fit,
    response.count = count.fit,
    par = par,
    theta = istem_random_theta(response.group$G, D, theta_mu, method$theta.lower,
                               method$theta.upper),
    Corr = istem_random_corr(D, include_corr),
    theta_mu = theta_mu,
    Q.matrix = Q.matrix,
    length.poly = length.poly,
    max.poly = max_poly,
    optim.maxit = method$optim.maxit,
    bounds = bounds,
    prior = prior
  )

  run <- istem_run(
    state = state,
    method = method,
    N = N,
    vis = vis,
    label = "MGPCM",
    update_parameters = mgpcm_istem_update_parameters,
    param_vec = function(state) {
      mgpcm_istem_param_vec(state$par, state$Corr, state$Q.matrix,
                            state$length.poly, include_corr)
    },
    sample_theta = mgpcm_istem_sample_theta,
    logLik_fun = function(state) {
      mgpcm_istem_loglik_trace(
        state, method$L, method$theta.lower, method$theta.upper
      )
    }
  )

  final <- mgpcm_istem_param_unpack(
    pv = run$chain.mean,
    pv_se = run$chain.sd,
    pv_rhat = NULL,
    I = I,
    D = D,
    Q.matrix = Q.matrix,
    length.poly = length.poly,
    max_poly = max_poly,
    Corr = run$state$Corr,
    include_corr = include_corr
  )

  theta.est <- istem_expand_group_matrix(
    run$theta.est, response.group$group, rownames(response))
  theta.se <- istem_expand_group_matrix(
    run$theta.se, response.group$group, rownames(response))
  theta_est <- istem_named_theta(theta.est, theta.se, response, D)
  par <- final$par
  par.se <- final$par.se
  par.Rhat <- final$par.rhat
  par.free <- mgpcm_par_free_mask(Q.matrix, length.poly, max_poly)
  colnames(par) <- colnames(par.se) <- colnames(par.Rhat) <- colnames(par.free) <-
    c(paste0("a", seq_len(D)), paste0("d", 0:(max_poly - 1L)))
  rownames(par) <- rownames(par.se) <- rownames(par.Rhat) <- rownames(par.free) <-
    colnames(response)
  par.se[!par.free] <- NA_real_
  par.Rhat[!par.free] <- NA_real_
  npar <- sum(par.free) + corr_df

  Corr.Rhat <- final$Corr.rhat
  dimnames(Corr.Rhat) <- list(paste0("Dim.", seq_len(D)), paste0("Dim.", seq_len(D)))

  results <- list(
    npar = npar,
    method = "iStEM",
    theta = theta_est,
    par = list(est = par, se = par.se, Rhat = par.Rhat, free = par.free),
    Corr = list(est = final$Corr, se = final$Corr.se, Rhat = Corr.Rhat),
    stan.obj = NULL,
    MCMC.obj = NULL,
    Q.matrix = Q.matrix,
    length.poly = length.poly,
    call = call,
    arguments = list(
      response = response,
      D = D,
      Q.matrix = Q.matrix,
      length.poly = length.poly,
      method = "iStEM",
      cores = common.method$cores,
      vis = vis,
      seed = common.method$seed,
      control.model = control.model,
      control.method = fit_effective_method_control(
        control.method, method = c(method, bounds), common = common.method
      )
    ),
    iStEM = utils::modifyList(run$iStEM, list(prior = prior))
  )

  class(results) <- "MGPCM"
  results$logLik <- logLik.MGPCM(results, L = method$L,
                                  theta.low = method$theta.lower,
                                  theta.up = method$theta.upper)
  results
}

mgpcm_istem_prior_control <- function(control.model) {
  use_prior <- get_ctrl("use.prior", TRUE, control.model)
  if (!is.logical(use_prior) || length(use_prior) != 1L || is.na(use_prior)) {
    stop("'use.prior' in 'control.model' must be TRUE or FALSE.",
         call. = FALSE)
  }
  list(
    use.prior = use_prior,
    a.mu = istem_scalar(get_ctrl("a.mu", 0.25, control.model), "a.mu"),
    a.sigma = istem_scalar(get_ctrl("a.sigma", 0.25, control.model),
                           "a.sigma", positive = TRUE),
    d.mu = istem_scalar(get_ctrl("d.mu", 0, control.model), "d.mu"),
    d.sigma = istem_scalar(get_ctrl("d.sigma", 1, control.model),
                           "d.sigma", positive = TRUE)
  )
}

mgpcm_istem_check_bounds <- function(bounds) {
  for (nm in names(bounds)) {
    bounds[[nm]] <- istem_scalar(bounds[[nm]], nm, arg = "control.method")
  }
  if (bounds$a.lower <= 0 || bounds$a.lower >= bounds$a.upper) {
    stop("'a.lower' must be positive and smaller than 'a.upper'.",
         call. = FALSE)
  }
  if (bounds$d.lower >= bounds$d.upper) {
    stop("'d.lower' must be smaller than 'd.upper'.", call. = FALSE)
  }
  invisible(NULL)
}

mgpcm_istem_initial_par <- function(response, Q.matrix, length.poly, bounds,
                                    prior) {
  I <- ncol(response)
  D <- ncol(Q.matrix)
  max_poly <- max(length.poly)
  par <- matrix(NA_real_, I, D + max_poly)
  a <- matrix(0, I, D)
  a[Q.matrix == 1] <- istem_rlnorm_bounded(
    sum(Q.matrix == 1), prior$a.mu, prior$a.sigma,
    bounds$a.lower, bounds$a.upper
  )
  par[, seq_len(D)] <- a
  for (i in seq_len(I)) {
    Ki <- length.poly[i]
    d_i <- rep(NA_real_, max_poly)
    d_i[1L] <- 0
    if (Ki > 1L) {
      d_i[2:Ki] <- istem_rnorm_bounded(
        Ki - 1L, prior$d.mu, prior$d.sigma,
        bounds$d.lower, bounds$d.upper
      )
    }
    par[i, D + seq_len(max_poly)] <- d_i
  }
  par
}

mgpcm_istem_sample_theta <- function(state, theta_grid_length,
                                     theta_lower, theta_upper) {
  chol <- istem_corr_chol(state$Corr)
  cpp_gibbs_mgpcm_theta(
    theta_ = state$theta,
    par_ = state$par,
    response_ = state$response,
    length_poly_ = as.integer(state$length.poly),
    theta_mu_ = state$theta_mu,
    chol_corr_ = chol$chol,
    log_diag_sum = chol$log_diag_sum,
    step = theta_grid_length,
    lower = theta_lower,
    upper = theta_upper
  )
}

mgpcm_istem_update_parameters <- function(state) {
  for (i in seq_len(ncol(state$response))) {
    state$par[i, ] <- mgpcm_istem_fit_item(
      y = state$response[, i],
      theta = state$theta,
      par_i = state$par[i, ],
      q_i = state$Q.matrix[i, ],
      Ki = state$length.poly[i],
      weight = istem_state_weight(state),
      optim_maxit = state$optim.maxit,
      bounds = state$bounds,
      prior = state$prior
    )
  }
  state
}

mgpcm_istem_fit_item <- function(y, theta, par_i, q_i, Ki, weight, optim_maxit,
                                 bounds, prior) {
  D <- ncol(theta)
  active <- q_i == 1
  x0 <- c(par_i[seq_len(D)][active], par_i[D + 2:Ki])
  lower <- c(rep(bounds$a.lower, sum(active)), rep(bounds$d.lower, Ki - 1L))
  upper <- c(rep(bounds$a.upper, sum(active)), rep(bounds$d.upper, Ki - 1L))

  obj <- function(x) {
    item <- mgpcm_istem_x_to_item(x, par_i, active, D, Ki)
    ll <- cpp_mgpcm_item_loglik_weighted(
      theta = theta,
      a = item$a,
      d = item$d,
      response = as.integer(y),
      weight = weight
    )
    lp <- mgpcm_istem_item_log_prior(item, active, prior)
    if (!is.finite(lp)) return(Inf)
    -ll - lp
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
  item <- mgpcm_istem_x_to_item(opt$par, par_i, active, D, Ki)
  out <- par_i
  out[seq_len(D)] <- item$a
  out[D + seq_len(Ki)] <- item$d
  out
}

mgpcm_istem_x_to_item <- function(x, par_i, active, D, Ki) {
  a <- rep(0, D)
  idx <- 1L
  n_active <- sum(active)
  if (n_active > 0L) {
    a[active] <- x[idx:(idx + n_active - 1L)]
    idx <- idx + n_active
  }
  d <- rep(NA_real_, Ki)
  d[1L] <- 0
  if (Ki > 1L) {
    d[2:Ki] <- x[idx:(idx + Ki - 2L)]
  }
  out <- par_i
  out[seq_len(D)] <- a
  out[D + seq_len(Ki)] <- d
  list(par = out, a = a, d = d)
}

mgpcm_istem_item_log_prior <- function(item, active, prior) {
  if (!isTRUE(prior$use.prior)) return(0)
  lp <- 0
  a_active <- item$a[active]
  if (length(a_active) > 0L) {
    if (any(a_active <= 0)) return(-Inf)
    lp <- lp + sum(stats::dlnorm(
      a_active, meanlog = prior$a.mu, sdlog = prior$a.sigma, log = TRUE
    ))
  }
  if (length(item$d) > 1L) {
    lp <- lp + sum(stats::dnorm(
      item$d[-1L], mean = prior$d.mu, sd = prior$d.sigma, log = TRUE
    ))
  }
  lp
}

mgpcm_istem_param_vec <- function(par, Corr, Q.matrix, length.poly,
                                  include_corr) {
  D <- ncol(Q.matrix)
  pv <- par[, seq_len(D), drop = FALSE][Q.matrix == 1]
  for (i in seq_along(length.poly)) {
    Ki <- length.poly[i]
    if (Ki > 1L) {
      pv <- c(pv, par[i, D + 2:Ki])
    }
  }
  c(pv, istem_param_corr_vec(Corr, include_corr))
}

mgpcm_istem_param_unpack <- function(pv, pv_se, pv_rhat = NULL, I, D, Q.matrix,
                                     length.poly, max_poly, Corr, include_corr) {
  idx <- 1L
  par <- matrix(NA_real_, I, D + max_poly)
  par.se <- matrix(NA_real_, I, D + max_poly)
  par.rhat <- matrix(NA_real_, I, D + max_poly)
  has_rhat <- !is.null(pv_rhat)
  a <- matrix(0, I, D)
  a.se <- matrix(NA_real_, I, D)
  a.rhat <- matrix(NA_real_, I, D)
  n_active <- sum(Q.matrix == 1)
  if (n_active > 0L) {
    a[Q.matrix == 1] <- pv[idx:(idx + n_active - 1L)]
    a.se[Q.matrix == 1] <- pv_se[idx:(idx + n_active - 1L)]
    if (has_rhat) a.rhat[Q.matrix == 1] <- pv_rhat[idx:(idx + n_active - 1L)]
    idx <- idx + n_active
  }
  par[, seq_len(D)] <- a
  par.se[, seq_len(D)] <- a.se
  par.rhat[, seq_len(D)] <- a.rhat
  for (i in seq_len(I)) {
    Ki <- length.poly[i]
    par[i, D + 1L] <- 0
    par.se[i, D + 1L] <- NA_real_
    if (Ki > 1L) {
      n_d <- Ki - 1L
      par[i, D + 2:Ki] <- pv[idx:(idx + n_d - 1L)]
      par.se[i, D + 2:Ki] <- pv_se[idx:(idx + n_d - 1L)]
      if (has_rhat) par.rhat[i, D + 2:Ki] <- pv_rhat[idx:(idx + n_d - 1L)]
      idx <- idx + n_d
    }
  }
  corr <- istem_unpack_corr(pv, pv_se, idx, D, Corr, include_corr)
  corr.rhat <- istem_unpack_corr_rhat(pv_rhat, idx, D, include_corr)
  list(par = par, par.se = par.se, par.rhat = par.rhat,
       Corr = corr$Corr, Corr.se = corr$Corr.se, Corr.rhat = corr.rhat)
}
