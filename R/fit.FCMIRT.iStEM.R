################################# FCMIRT iStEM #################################

fit.FCMIRT.iStEM <- function(data, model = "m2pl", Q.matrix = NULL,
                             block.items = NULL, D = NULL, fc.type = "RANK",
                             control.model = NULL, control.method = NULL,
                             .call = NULL) {

  call <- if (is.null(.call)) match.call() else .call
  control.model <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  vis <- common.method$vis
  method <- istem_method_control(control.method)
  set.seed(common.method$seed)

  model <- tolower(model)
  model_type <- model_type_to_int(model)

  fc <- istem_prepare_fc_data(data, block.items, fc.type)
  I <- fc$I
  N <- fc$N
  response.group <- istem_response_groups(fc$response)
  response.fit <- response.group$response
  count.fit <- response.group$count
  if (!is.null(Q.matrix)) {
    Q.matrix <- as.matrix(Q.matrix)
    if (anyNA(Q.matrix) || !all(Q.matrix %in% c(0, 1))) {
      stop("'Q.matrix' must be a matrix containing only 0 or 1.",
           call. = FALSE)
    }
    if (nrow(Q.matrix) != I) {
      stop("'Q.matrix' must have ", I, " rows, one for each item.",
           call. = FALSE)
    }
    if (any(rowSums(Q.matrix) < 1L)) {
      stop("Each item must measure at least one trait (each row of 'Q.matrix' needs a 1).",
           call. = FALSE)
    }
    D.q <- ncol(Q.matrix)
    if (!is.null(D) && D != D.q) {
      warning("'D' (", D, ") does not match ncol(Q.matrix) (", D.q,
              "); using ncol(Q.matrix).")
    }
    D <- D.q
  }
  if (is.null(D) || is.na(D)) {
    stop("Cannot determine 'D'. Please provide either 'Q.matrix' or 'D'.",
         call. = FALSE)
  }
  D <- as.integer(D)
  if (D < 1L) {
    stop("'D' must be a positive integer.", call. = FALSE)
  }
  if (is.null(Q.matrix)) {
    Q.matrix <- matrix(1, I, D)
  }
  storage.mode(Q.matrix) <- "numeric"
  method <- istem_apply_grid_control(method, control.model, D)
  istem_check_method_control(method)
  control.model <- istem_effective_grid_model_control(control.model, method)

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

  theta_mu <- array(get_ctrl("theta.mu", rep(0, D), control.model), dim = D)
  include_corr <- D > 1L && !method$fix.corr
  par <- fcmirt_istem_initial_par(model_type, Q.matrix, bounds, prior,
                                  fc$block.items)

  state <- list(
    response = response.fit,
    response.count = count.fit,
    par = par,
    theta = istem_random_theta(response.group$G, D, theta_mu, method$theta.lower,
                               method$theta.upper),
    Corr = istem_random_corr(D, include_corr),
    theta_mu = theta_mu,
    model.type = model_type,
    Q.matrix = Q.matrix,
    block.items = fc$block.items,
    patterns = fc$patterns,
    patterns.total = fc$patterns.total,
    optim.maxit = method$optim.maxit,
    bounds = bounds,
    prior = prior
  )

  run <- istem_run(
    state = state,
    method = method,
    N = N,
    vis = vis,
    label = "FCMIRT",
    update_parameters = fcmirt_istem_update_parameters,
    param_vec = function(state) {
      fcmirt_istem_param_vec(state$par, state$Corr, state$model.type,
                             state$Q.matrix, state$block.items, include_corr)
    },
    sample_theta = fcmirt_istem_sample_theta,
    logLik_fun = function(state) {
      fcmirt_istem_loglik_trace(
        state, method$L, method$theta.lower, method$theta.upper
      )
    }
  )

  final <- fcmirt_istem_param_unpack(
    pv = run$chain.mean,
    pv_se = run$chain.sd,
    pv_rhat = NULL,
    I = I, D = D,
    model_type = model_type,
    Q.matrix = Q.matrix,
    block.items = fc$block.items,
    Corr = run$state$Corr,
    include_corr = include_corr
  )

  theta.est <- istem_expand_group_matrix(
    run$theta.est, response.group$group, rownames(fc$response))
  theta.se <- istem_expand_group_matrix(
    run$theta.se, response.group$group, rownames(fc$response))
  theta_est <- istem_named_theta(theta.est, theta.se, fc$response, D)
  par <- final$par
  par.se <- final$par.se
  par.Rhat <- final$par.rhat
  par.free <- fcmirt_par_free_mask(model_type, Q.matrix, fc$block.items)
  colnames(par) <- colnames(par.se) <- colnames(par.Rhat) <- colnames(par.free) <-
    c(paste0("a", seq_len(D)), "b", "c", "d")
  rownames(par) <- rownames(par.se) <- rownames(par.Rhat) <- rownames(par.free) <-
    paste0("item", fc$all.items)
  par.se[!par.free] <- NA_real_
  par.Rhat[!par.free] <- NA_real_

  corr_df <- free_corr_npar(D, include_corr)
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
    model = model,
    Q.matrix = Q.matrix,
    block.items = fc$block.items,
    fc.type = fc$fc.type,
    response = fc$response,
    patterns = fc$patterns,
    patterns.total = fc$patterns.total,
    call = call,
    arguments = list(
      data = data,
      model = model,
      Q.matrix = Q.matrix,
      block.items = fc$block.items,
      D = D,
      fc.type = fc$fc.type,
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

  class(results) <- "FCMIRT"
  results$logLik <- logLik.FCMIRT(results, L = method$L,
                                    theta.low = method$theta.lower,
                                    theta.up = method$theta.upper)
  results
}

fcmirt_par_free_mask <- function(model_type, Q.matrix, block.items) {
  I <- nrow(Q.matrix)
  D <- ncol(Q.matrix)
  free <- matrix(FALSE, nrow = I, ncol = D + 3L)
  if (model_type >= 2L) {
    free[, seq_len(D)] <- Q.matrix == 1
  }
  for (items in block.items) {
    K <- length(items)
    free[items[-K], D + 1L] <- TRUE
  }
  if (model_type >= 3L) {
    free[, D + 2L] <- TRUE
  }
  if (model_type == 4L) {
    free[, D + 3L] <- TRUE
  }
  colnames(free) <- c(paste0("a", seq_len(D)), "b", "c", "d")
  free
}

fcmirt_istem_initial_par <- function(model_type, Q.matrix, bounds, prior,
                                     block.items) {
  I <- nrow(Q.matrix)
  D <- ncol(Q.matrix)
  par <- matrix(NA_real_, I, D + 3L)
  if (model_type == 1L) {
    par[, seq_len(D)] <- Q.matrix
  } else {
    a <- matrix(0, I, D)
    a[Q.matrix == 1] <- istem_rlnorm_bounded(
      sum(Q.matrix == 1), prior$a.mu, prior$a.sigma,
      bounds$a.lower, bounds$a.upper
    )
    par[, seq_len(D)] <- a
  }
  par[, D + 1L] <- NA_real_
  for (items in block.items) {
    b <- istem_rnorm_bounded(length(items), prior$b.mu, prior$b.sigma,
                             bounds$b.lower, bounds$b.upper)
    b <- istem_clip(b - mean(b), bounds$b.lower, bounds$b.upper)
    b <- b - mean(b)
    par[items, D + 1L] <- b
  }
  par[, D + 2L] <- if (model_type >= 3L) {
    stats::runif(I, max(bounds$c.lower, prior$c.lower),
                 min(bounds$c.upper, prior$c.upper))
  } else {
    0
  }
  par[, D + 3L] <- if (model_type == 4L) {
    stats::runif(I, max(bounds$d.lower, prior$d.lower),
                 min(bounds$d.upper, prior$d.upper))
  } else {
    1
  }
  par
}

fcmirt_istem_sample_theta <- function(state, theta_grid_length,
                                      theta_lower, theta_upper) {
  chol <- istem_corr_chol(state$Corr)
  cpp_gibbs_fcmirt_theta(
    theta_ = state$theta,
    par_ = state$par,
    response_ = state$response,
    theta_mu_ = state$theta_mu,
    chol_corr_ = chol$chol,
    log_diag_sum = chol$log_diag_sum,
    patterns_total = state$patterns.total,
    patterns = state$patterns,
    step = theta_grid_length,
    lower = theta_lower,
    upper = theta_upper
  )
}

fcmirt_istem_update_parameters <- function(state) {
  for (b in seq_along(state$block.items)) {
    items <- state$block.items[[b]]
    state$par[items, ] <- fcmirt_istem_fit_block(
      y = state$response[, b],
      theta = state$theta,
      par_block = state$par[items, , drop = FALSE],
      q_block = state$Q.matrix[items, , drop = FALSE],
      patterns_total = cpp_istem_local_patterns(
        as.matrix(state$patterns.total[[b]]), as.integer(items)
      ),
      patterns = cpp_istem_local_patterns(
        as.matrix(state$patterns[[b]]), as.integer(items)
      ),
      model_type = state$model.type,
      weight = istem_state_weight(state),
      optim_maxit = state$optim.maxit,
      bounds = state$bounds,
      prior = state$prior
    )
  }
  state
}

fcmirt_istem_fit_block <- function(y, theta, par_block, q_block,
                                   patterns_total, patterns, model_type,
                                   weight, optim_maxit, bounds, prior) {
  D <- ncol(q_block)
  active <- q_block == 1
  x0 <- numeric(0)
  lower <- numeric(0)
  upper <- numeric(0)
  if (model_type >= 2L) {
    x0 <- c(x0, par_block[, seq_len(D), drop = FALSE][active])
    lower <- c(lower, rep(bounds$a.lower, sum(active)))
    upper <- c(upper, rep(bounds$a.upper, sum(active)))
  }
  x0 <- c(x0, par_block[-nrow(par_block), D + 1L])
  lower <- c(lower, rep(bounds$b.lower, nrow(par_block) - 1L))
  upper <- c(upper, rep(bounds$b.upper, nrow(par_block) - 1L))
  if (model_type >= 3L) {
    x0 <- c(x0, par_block[, D + 2L])
    lower <- c(lower, rep(bounds$c.lower, nrow(par_block)))
    upper <- c(upper, rep(bounds$c.upper, nrow(par_block)))
  }
  if (model_type == 4L) {
    x0 <- c(x0, par_block[, D + 3L])
    lower <- c(lower, rep(bounds$d.lower, nrow(par_block)))
    upper <- c(upper, rep(bounds$d.upper, nrow(par_block)))
  }

  obj <- function(x) {
    block <- fcmirt_istem_x_to_block(x, par_block, q_block, model_type)
    K <- nrow(block)
    if (K > 1L && any(block[-K, D + 1L] < bounds$b.lower |
                       block[-K, D + 1L] > bounds$b.upper)) {
      return(Inf)
    }
    b_last <- block[K, D + 1L]
    penalty <- 0
    if (b_last < bounds$b.lower || b_last > bounds$b.upper) {
      excess <- max(bounds$b.lower - b_last, b_last - bounds$b.upper, 0)
      penalty <- 10 * excess^2
    }
    ll <- cpp_fcmirt_block_loglik_weighted(
      theta = theta,
      par_block = block,
      response = as.integer(y),
      patterns_total = patterns_total,
      patterns = patterns,
      weight = weight
    )
    lp <- fcmirt_istem_block_log_prior(block, q_block, model_type, prior)
    if (!is.finite(lp)) return(Inf)
    -ll - lp + penalty
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
    return(par_block)
  }
  fcmirt_istem_x_to_block(opt$par, par_block, q_block, model_type)
}

fcmirt_istem_x_to_block <- function(x, par_block, q_block, model_type) {
  D <- ncol(q_block)
  I <- nrow(q_block)
  idx <- 1L
  out <- par_block
  if (model_type == 1L) {
    out[, seq_len(D)] <- q_block
  } else {
    a <- matrix(0, I, D)
    n_active <- sum(q_block == 1)
    if (n_active > 0L) {
      a[q_block == 1] <- x[idx:(idx + n_active - 1L)]
      idx <- idx + n_active
    }
    out[, seq_len(D)] <- a
  }
  b_free <- x[idx:(idx + I - 2L)]
  out[-I, D + 1L] <- b_free
  out[I, D + 1L] <- -sum(b_free)
  idx <- idx + I - 1L
  out[, D + 2L] <- if (model_type >= 3L) {
    val <- x[idx:(idx + I - 1L)]
    idx <- idx + I
    val
  } else {
    0
  }
  out[, D + 3L] <- if (model_type == 4L) {
    x[idx:(idx + I - 1L)]
  } else {
    1
  }
  out
}

fcmirt_istem_block_log_prior <- function(par_block, q_block, model_type, prior) {
  if (!isTRUE(prior$use.prior)) return(0)
  D <- ncol(q_block)
  lp <- sum(stats::dnorm(
    par_block[, D + 1L], mean = prior$b.mu, sd = prior$b.sigma, log = TRUE
  ))
  if (model_type >= 2L) {
    a_active <- par_block[, seq_len(D), drop = FALSE][q_block == 1]
    if (any(a_active <= 0)) return(-Inf)
    lp <- lp + sum(stats::dlnorm(
      a_active, meanlog = prior$a.mu, sdlog = prior$a.sigma, log = TRUE
    ))
  }
  if (model_type >= 3L &&
      any(par_block[, D + 2L] < prior$c.lower |
          par_block[, D + 2L] > prior$c.upper)) {
    return(-Inf)
  }
  if (model_type == 4L &&
      any(par_block[, D + 3L] < prior$d.lower |
          par_block[, D + 3L] > prior$d.upper)) {
    return(-Inf)
  }
  lp
}

fcmirt_istem_param_vec <- function(par, Corr, model_type, Q.matrix,
                                   block.items, include_corr) {
  D <- ncol(Q.matrix)
  b_free <- unlist(lapply(block.items, function(items) {
    par[items[-length(items)], D + 1L]
  }), use.names = FALSE)
  pv <- if (model_type == 1L) {
    b_free
  } else {
    c(par[, seq_len(D), drop = FALSE][Q.matrix == 1], b_free)
  }
  if (model_type >= 3L) {
    pv <- c(pv, par[, D + 2L])
  }
  if (model_type == 4L) {
    pv <- c(pv, par[, D + 3L])
  }
  c(pv, istem_param_corr_vec(Corr, include_corr))
}

fcmirt_istem_param_unpack <- function(pv, pv_se, pv_rhat = NULL, I, D, model_type,
                                      Q.matrix, block.items, Corr, include_corr) {
  idx <- 1L
  par <- matrix(NA_real_, I, D + 3L)
  par.se <- matrix(NA_real_, I, D + 3L)
  par.rhat <- matrix(NA_real_, I, D + 3L)
  has_rhat <- !is.null(pv_rhat)
  if (model_type == 1L) {
    par[, seq_len(D)] <- Q.matrix
    par.se[, seq_len(D)] <- NA_real_
  } else {
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
  }
  par[, D + 1L] <- 0
  par.se[, D + 1L] <- NA_real_
  for (items in block.items) {
    K <- length(items)
    b_free <- pv[idx:(idx + K - 2L)]
    b_free_se <- pv_se[idx:(idx + K - 2L)]
    par[items[-K], D + 1L] <- b_free
    par[items[K], D + 1L] <- -sum(b_free)
    par.se[items[-K], D + 1L] <- b_free_se
    if (has_rhat) {
      b_free_rhat <- pv_rhat[idx:(idx + K - 2L)]
      par.rhat[items[-K], D + 1L] <- b_free_rhat
    }
    idx <- idx + K - 1L
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
  corr <- istem_unpack_corr(pv, pv_se, idx, D, Corr, include_corr)
  corr.rhat <- if (has_rhat) {
    istem_unpack_corr_rhat(pv_rhat, idx, D, include_corr)
  } else {
    matrix(NA_real_, D, D)
  }
  list(par = par, par.se = par.se, par.rhat = par.rhat,
       Corr = corr$Corr, Corr.se = corr$Corr.se, Corr.rhat = corr.rhat)
}
