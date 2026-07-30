################################# FCGGUM iStEM #################################

fit.FCGGUM.iStEM <- function(data, Q.matrix = NULL, block.items = NULL,
                             D = NULL, fc.type = "RANK",
                             control.model = NULL, control.method = NULL,
                             .call = NULL) {

  call <- if (is.null(.call)) match.call() else .call
  control.model <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  vis <- common.method$vis
  method <- istem_method_control(control.method)
  set.seed(common.method$seed)

  fc <- istem_prepare_fc_data(data, block.items, fc.type)
  I <- fc$I
  N <- fc$N
  response.group <- istem_response_groups(fc$response)
  response.fit <- response.group$response
  count.fit <- response.group$count
  if (!is.null(Q.matrix)) {
    Q.matrix <- as.matrix(Q.matrix)
    if (anyNA(Q.matrix) || !all(Q.matrix %in% c(-1, 0, 1))) {
      stop("'Q.matrix' must be a matrix containing only -1, 0, or 1.",
           call. = FALSE)
    }
    if (nrow(Q.matrix) != I) {
      stop("'Q.matrix' must have ", I, " rows, one for each item.",
           call. = FALSE)
    }
    if (any(rowSums(abs(Q.matrix)) < 1L)) {
      stop("Each item must measure at least one trait (each row of 'Q.matrix' needs a non-zero entry).",
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
  if (D < 2L) {
    stop("FCGGUM requires at least 2 latent dimensions (D >= 2).",
         call. = FALSE)
  }
  if (is.null(Q.matrix)) {
    Q.matrix <- matrix(1, I, D)
  }
  storage.mode(Q.matrix) <- "numeric"
  method <- istem_apply_grid_control(method, control.model, D)
  istem_check_method_control(method)
  control.model <- istem_effective_grid_model_control(control.model, method)

  length.poly <- rep(2L, I)
  max_poly <- 2L
  prior <- mggum_istem_prior_control(control.model)
  bounds <- list(
    a.lower = get_ctrl("a.lower", 1e-4, control.method),
    a.upper = get_ctrl("a.upper", 6, control.method),
    delta.lower = get_ctrl("delta.lower", -6, control.method),
    delta.upper = get_ctrl("delta.upper", 6, control.method),
    tau.gap.lower = get_ctrl("tau.gap.lower", 1e-4, control.method),
    tau.gap.upper = get_ctrl("tau.gap.upper", 6, control.method)
  )
  mggum_istem_check_bounds(bounds)

  theta_mu <- array(get_ctrl("theta.mu", rep(0, D), control.model), dim = D)
  include_corr <- D > 1L && !method$fix.corr
  corr_df <- free_corr_npar(D, include_corr)
  par <- mggum_istem_initial_par(Q.matrix, length.poly, max_poly, bounds,
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
    label = "FCGGUM",
    update_parameters = fcggum_istem_update_parameters,
    param_vec = function(state) {
      mggum_istem_param_vec(state$par, state$Corr, state$Q.matrix,
                            state$length.poly, include_corr)
    },
    sample_theta = fcggum_istem_sample_theta,
    logLik_fun = function(state) {
      fcggum_istem_loglik_trace(
        state, method$L, method$theta.lower, method$theta.upper
      )
    }
  )

  final <- mggum_istem_param_unpack(
    pv = run$chain.mean,
    pv_se = run$chain.sd,
    pv_rhat = NULL,
    I = I, D = D,
    Q.matrix = Q.matrix,
    length.poly = length.poly,
    max_poly = max_poly,
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
  par.free <- mggum_par_free_mask(Q.matrix, length.poly, max_poly)
  colnames(par) <- colnames(par.se) <- colnames(par.Rhat) <- colnames(par.free) <-
    c(paste0("a", seq_len(D)),
      paste0("delta", seq_len(D)),
      paste0("tau", 0:(max_poly - 1L)))
  rownames(par) <- rownames(par.se) <- rownames(par.Rhat) <- rownames(par.free) <-
    paste0("item", fc$all.items)
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
    block.items = fc$block.items,
    fc.type = fc$fc.type,
    response = fc$response,
    length.poly = length.poly,
    patterns = fc$patterns,
    patterns.total = fc$patterns.total,
    call = call,
    arguments = list(
      data = data,
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

  class(results) <- "FCGGUM"
  results$logLik <- logLik.FCGGUM(results, L = method$L,
                                    theta.low = method$theta.lower,
                                    theta.up = method$theta.upper)
  results
}

fcggum_istem_sample_theta <- function(state, theta_grid_length,
                                      theta_lower, theta_upper) {
  chol <- istem_corr_chol(state$Corr)
  cpp_gibbs_fcggum_theta(
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

fcggum_istem_update_parameters <- function(state) {
  for (b in seq_along(state$block.items)) {
    items <- state$block.items[[b]]
    state$par[items, ] <- fcggum_istem_fit_block(
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
      weight = istem_state_weight(state),
      optim_maxit = state$optim.maxit,
      bounds = state$bounds,
      prior = state$prior
    )
  }
  state
}

fcggum_istem_fit_block <- function(y, theta, par_block, q_block,
                                   patterns_total, patterns, weight, optim_maxit,
                                   bounds, prior) {
  D <- ncol(q_block)
  active <- q_block != 0
  x0 <- c(par_block[, seq_len(D), drop = FALSE][active],
          par_block[, D + seq_len(D), drop = FALSE][active],
          -par_block[, D + D + 2L])
  delta_lower <- ifelse(q_block[active] > 0, max(0, bounds$delta.lower),
                        bounds$delta.lower)
  delta_upper <- ifelse(q_block[active] > 0, bounds$delta.upper,
                        min(0, bounds$delta.upper))
  lower <- c(rep(bounds$a.lower, sum(active)),
             delta_lower,
             rep(bounds$tau.gap.lower, nrow(par_block)))
  upper <- c(rep(bounds$a.upper, sum(active)),
             delta_upper,
             rep(bounds$tau.gap.upper, nrow(par_block)))

  obj <- function(x) {
    block <- fcggum_istem_x_to_block(x, par_block, q_block)
    ll <- cpp_fcggum_block_loglik_weighted(
      theta = theta,
      par_block = block,
      response = as.integer(y),
      patterns_total = patterns_total,
      patterns = patterns,
      weight = weight
    )
    lp <- 0
    for (i in seq_len(nrow(block))) {
      item <- list(
        a = block[i, seq_len(D)],
        delta = block[i, D + seq_len(D)],
        tau = block[i, D + D + 1:2]
      )
      lp <- lp + mggum_istem_item_log_prior(item, q_block[i, ], prior)
    }
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
    return(par_block)
  }
  fcggum_istem_x_to_block(opt$par, par_block, q_block)
}

fcggum_istem_x_to_block <- function(x, par_block, q_block) {
  D <- ncol(q_block)
  I <- nrow(q_block)
  active <- q_block != 0
  idx <- 1L
  out <- par_block
  a <- matrix(0, I, D)
  n_active <- sum(active)
  if (n_active > 0L) {
    a[active] <- x[idx:(idx + n_active - 1L)]
    idx <- idx + n_active
  }
  delta <- matrix(0, I, D)
  if (n_active > 0L) {
    delta[active] <- x[idx:(idx + n_active - 1L)]
    idx <- idx + n_active
  }
  tau_gap <- x[idx:(idx + I - 1L)]
  out[, seq_len(D)] <- a
  out[, D + seq_len(D)] <- delta
  out[, D + D + 1L] <- 0
  out[, D + D + 2L] <- -tau_gap
  out
}
