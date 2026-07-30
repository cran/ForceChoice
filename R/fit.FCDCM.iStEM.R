fit.FCDCM.iStEM <- function(data, Q.matrix, block.items = NULL,
                            dcm.type = "DINA",
                            control.model = NULL,
                            control.method = NULL,
                            .call = NULL) {
  call <- if (is.null(.call)) match.call() else .call
  control.model <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  method <- istem_method_control(control.method)
  set.seed(common.method$seed)

  data.matrix <- as.matrix(data)
  N <- nrow(data.matrix)
  B <- ncol(data.matrix)
  if (is.null(N) || is.null(B) || N < 1L || B < 1L) {
    stop("'data' must be an N x B FCDCM data matrix.", call. = FALSE)
  }
  if (is.numeric(data.matrix) || is.integer(data.matrix) ||
      is.logical(data.matrix)) {
    stop("'data' must contain FCDCM ranking strings such as '1>3'. ",
         "Use get.data.from.response.FCDCM() to convert a response matrix.",
         call. = FALSE)
  }

  if (is.null(block.items)) {
    block.items <- get.block.items.from.data.FCDCM(data)
  }
  if (!is.list(block.items) || length(block.items) != B) {
    stop("'block.items' must be a list with one two-statement entry per block.",
         call. = FALSE)
  }
  block.items <- lapply(block.items, as.integer)
  if (any(vapply(block.items, length, integer(1L)) != 2L)) {
    stop("Each FCDCM block must contain exactly two statement indices.",
         call. = FALSE)
  }
  all.items <- unlist(block.items, use.names = FALSE)
  if (anyNA(all.items) || any(all.items < 1L)) {
    stop("FCDCM statement indices in 'block.items' must be positive integers.",
         call. = FALSE)
  }
  if (anyDuplicated(all.items)) {
    stop("'block.items' contains duplicate statement indices.", call. = FALSE)
  }
  I.states <- length(all.items)
  if (!setequal(all.items, seq_len(I.states))) {
    stop("'block.items' must contain exactly the statement indices 1 through ",
         I.states, ".", call. = FALSE)
  }

  patterns <- matrix(all.items, ncol = 2L, byrow = TRUE)
  storage.mode(patterns) <- "integer"

  response <- get.response.from.data.FCDCM(data.matrix, block.items)$response
  if (anyNA(response) || !all(response %in% c(0L, 1L))) {
    stop("FCDCM responses must contain only 0/1 after parsing.",
         call. = FALSE)
  }
  storage.mode(response) <- "integer"
  response.group <- istem_response_groups(response)
  response.fit <- response.group$response
  count.fit <- response.group$count

  Q.matrix <- as.matrix(Q.matrix)
  if (anyNA(Q.matrix) || !all(Q.matrix %in% c(0, 1))) {
    stop("'Q.matrix' must be a 0/1 matrix.", call. = FALSE)
  }
  if (nrow(Q.matrix) != I.states) {
    stop("'Q.matrix' must have ", I.states,
         " rows, one for each statement.", call. = FALSE)
  }
  if (any(rowSums(Q.matrix) < 1L)) {
    stop("Each statement must require at least one attribute.",
         call. = FALSE)
  }
  D <- ncol(Q.matrix)
  storage.mode(Q.matrix) <- "numeric"
  method <- istem_apply_grid_control(method, control.model, 1L)
  istem_check_method_control(method)
  control.model <- istem_effective_grid_model_control(control.model, method)

  dcm.type <- toupper(dcm.type)
  if (length(dcm.type) != I.states) {
    dcm.type <- rep(dcm.type[1L], I.states)
  }
  if (!all(dcm.type %in% c("DINA", "DINO"))) {
    stop("'dcm.type' must contain only 'DINA' or 'DINO'.",
         call. = FALSE)
  }

  alpha.patterns <- as.matrix(do.call(expand.grid, rep(list(0:1), D)))
  storage.mode(alpha.patterns) <- "numeric"
  zeta.patterns <- get.zeta(alpha.patterns, Q.matrix, dcm.type)
  storage.mode(zeta.patterns) <- "integer"

  block.names <- colnames(response)
  if (is.null(block.names)) block.names <- paste0("B", seq_len(B))
  dim.names <- colnames(Q.matrix)
  if (is.null(dim.names)) dim.names <- paste0("Dim.", seq_len(D))
  item.names <- rownames(Q.matrix)
  if (is.null(item.names)) item.names <- paste0("S", seq_len(I.states))

  colnames(response) <- block.names
  rownames(patterns) <- block.names
  colnames(patterns) <- c("SA", "SB")
  rownames(Q.matrix) <- item.names
  colnames(Q.matrix) <- dim.names
  colnames(alpha.patterns) <- dim.names
  colnames(zeta.patterns) <- item.names

  prep <- list(
    response = response,
    response.unique = response.fit,
    response.count = count.fit,
    response.group = response.group$group,
    G = response.group$G,
    Q.matrix = Q.matrix,
    block.items = block.items,
    patterns = patterns,
    dcm.type = dcm.type,
    alpha.patterns = alpha.patterns,
    zeta.patterns = zeta.patterns,
    N = N,
    B = B,
    D = D,
    I.states = I.states,
    block.names = block.names,
    dim.names = dim.names,
    item.names = item.names
  )
  par.init <- get_ctrl("par", NULL, control.model)

  prior <- list(
    use.prior = get_ctrl("use.prior", TRUE, control.model),
    delta1.mu = istem_scalar(get_ctrl("delta1.mu", 0.0, control.model),
                             "delta1.mu"),
    delta1.sigma = istem_scalar(get_ctrl("delta1.sigma", 0.50, control.model),
                                "delta1.sigma", positive = TRUE),
    delta0.mu = istem_scalar(get_ctrl("delta0.mu", 0.00, control.model),
                             "delta0.mu"),
    delta0.sigma = istem_scalar(get_ctrl("delta0.sigma", 1.00, control.model),
                                "delta0.sigma", positive = TRUE),
    eta0.mu = istem_scalar(
      get_ctrl("eta0.mu", 0.10, control.model),
      "eta0.mu", positive = TRUE),
    eta0.sigma = istem_scalar(
      get_ctrl("eta0.sigma", 0.10, control.model),
      "eta0.sigma", positive = TRUE),
    etaAB.mu = istem_scalar(
      get_ctrl("etaAB.mu", 0.30, control.model),
      "etaAB.mu", positive = TRUE),
    etaAB.sigma = istem_scalar(
      get_ctrl("etaAB.sigma", 0.10, control.model),
      "etaAB.sigma", positive = TRUE)
  )
  if (!is.logical(prior$use.prior) || length(prior$use.prior) != 1L ||
      is.na(prior$use.prior)) {
    stop("'use.prior' in 'control.model' must be TRUE or FALSE.",
         call. = FALSE)
  }

  if (is.null(par.init)) {
    par.init <- cbind(
      eta0  = rep(prior$eta0.mu, prep$B),
      etaAB = rep(prior$etaAB.mu, prep$B)
    )
  }
  par.init <- fcdcm_normalize_par(par.init, B = prep$B, arg = "par")
  rownames(par.init) <- prep$block.names

  delta1.init <- as.numeric(get_ctrl(
    "delta1",
    stats::rlnorm(prep$D, meanlog = prior$delta1.mu,
                  sdlog = prior$delta1.sigma),
    control.model
  ))
  delta0.init <- as.numeric(get_ctrl(
    "delta0",
    stats::rnorm(prep$D, mean = prior$delta0.mu, sd = prior$delta0.sigma),
    control.model
  ))
  if (length(delta1.init) != prep$D || length(delta0.init) != prep$D ||
      anyNA(delta1.init) || anyNA(delta0.init) ||
      any(!is.finite(delta1.init)) || any(!is.finite(delta0.init)) ||
      any(delta1.init <= 0)) {
    stop("'delta1' and 'delta0' initial values must be finite vectors of length D, with delta1 > 0.",
         call. = FALSE)
  }
  names(delta1.init) <- names(delta0.init) <- prep$dim.names

  theta.init <- get_ctrl(
    "theta",
    drop(istem_random_theta(prep$G, 1L, 0, method$theta.lower,
                            method$theta.upper)),
    control.model
  )
  theta.init <- as.numeric(theta.init)
  if (length(theta.init) != prep$G || anyNA(theta.init) ||
      any(!is.finite(theta.init))) {
    stop("'theta' initial values must be a finite vector of length G.",
         call. = FALSE)
  }

  class.prob.init <- model.FCDCM(prep$zeta.patterns, par.init,
                                 prep$patterns)
  alpha.init <- cpp_fcdcm_sample_alpha(
    theta = matrix(theta.init, ncol = 1L),
    response = prep$response.unique,
    delta1 = delta1.init,
    delta0 = delta0.init,
    alpha_patterns = prep$alpha.patterns,
    zeta_patterns = prep$zeta.patterns,
    class_prob = class.prob.init
  )

  state <- list(
    response = prep$response.unique,
    response.count = prep$response.count,
    theta = matrix(theta.init, ncol = 1L),
    alpha = alpha.init$alpha,
    zeta = alpha.init$zeta,
    class.index = alpha.init$class.index,
    Corr = matrix(1, 1L, 1L),
    delta1 = delta1.init,
    delta0 = delta0.init,
    par = par.init,
    alpha.patterns = prep$alpha.patterns,
    zeta.patterns = prep$zeta.patterns,
    patterns = prep$patterns,
    class.prob = class.prob.init,
    optim.maxit = method$optim.maxit,
    prior = prior
  )

  run <- istem_run(
    state = state,
    method = method,
    N = prep$N,
    vis = common.method$vis,
    label = "FCDCM",
    update_parameters = fcdcm_istem_update_parameters,
    param_vec = fcdcm_istem_param_vec,
    sample_theta = fcdcm_istem_sample_theta,
    logLik_fun = function(state) {
      fcdcm_istem_loglik_trace(
        state, method$L, method$theta.lower, method$theta.upper
      )
    }
  )

  final <- fcdcm_istem_param_unpack(
    pv = run$chain.mean,
    pv_se = run$chain.sd,
    pv_rhat = NULL,
    D = prep$D,
    B = prep$B,
    dim.names = prep$dim.names,
    block.names = prep$block.names,
    par.names = colnames(par.init)
  )

  posterior <- fcdcm_marginal_posterior(
    response = prep$response,
    delta1 = final$delta[, "delta1"],
    delta0 = final$delta[, "delta0"],
    par = final$par,
    alpha.patterns = prep$alpha.patterns,
    zeta.patterns = prep$zeta.patterns,
    patterns = prep$patterns,
    L = method$L,
    theta.lower = method$theta.lower,
    theta.upper = method$theta.upper
  )
  class.post <- posterior$class.post
  theta <- posterior$theta
  colnames(theta) <- "theta"
  rownames(theta) <- rownames(prep$response)
  theta.se <- posterior$theta.se
  colnames(theta.se) <- "theta"
  rownames(theta.se) <- rownames(prep$response)
  theta.Rhat <- matrix(NA_real_, prep$N, 1L, dimnames = dimnames(theta))

  alpha.prob <- class.post %*% prep$alpha.patterns
  colnames(alpha.prob) <- prep$dim.names
  rownames(alpha.prob) <- rownames(prep$response)
  alpha.est <- fcdcm_alpha_map(class.post, prep$alpha.patterns)
  storage.mode(alpha.est) <- "integer"
  alpha.se <- matrix(NA_real_, prep$N, prep$D,
                     dimnames = dimnames(alpha.est))
  alpha.Rhat <- alpha.se

  par.free <- matrix(TRUE, prep$B, ncol(final$par),
                     dimnames = list(prep$block.names,
                                     colnames(final$par)))
  npar <- ncol(final$par) * prep$B + 2L * prep$D

  results <- list(
    npar = npar,
    method = "iStEM",
    theta = list(est = theta, se = theta.se, Rhat = theta.Rhat),
    alpha = list(est = alpha.est, se = alpha.se, Rhat = alpha.Rhat,
                 prob = alpha.prob),
    class.post = class.post,
    delta = list(est = final$delta, se = final$delta.se,
                 Rhat = final$delta.Rhat),
    par = list(est = final$par, se = final$par.se,
               Rhat = final$par.Rhat, free = par.free),
    Corr = list(
      est = matrix(1, 1L, 1L, dimnames = list("theta", "theta")),
      se = matrix(0, 1L, 1L, dimnames = list("theta", "theta")),
      Rhat = matrix(NA_real_, 1L, 1L, dimnames = list("theta", "theta"))
    ),
    stan.obj = NULL,
    MCMC.obj = NULL,
    Q.matrix = prep$Q.matrix,
    block.items = prep$block.items,
    patterns = prep$patterns,
    alpha.patterns = prep$alpha.patterns,
    zeta.patterns = prep$zeta.patterns,
    dcm.type = prep$dcm.type,
    response = prep$response,
    call = call,
    arguments = list(
      data = data,
      Q.matrix = prep$Q.matrix,
      block.items = prep$block.items,
      D = prep$D,
      dcm.type = prep$dcm.type,
      method = "iStEM",
      cores = common.method$cores,
      vis = common.method$vis,
      seed = common.method$seed,
      control.model = control.model,
      control.method = fit_effective_method_control(
        control.method, method = method, common = common.method
      )
    ),
    iStEM = utils::modifyList(run$iStEM, list(prior = prior))
  )

  class(results) <- "FCDCM"
  results$logLik <- logLik.FCDCM(results, L = method$L,
                                  theta.low = method$theta.lower,
                                  theta.up = method$theta.upper)
  results
}

fcdcm_istem_sample_theta <- function(state, theta_grid_length,
                                     theta_lower, theta_upper) {
  cpp_gibbs_fcdcm_theta(
    theta = state$theta,
    response = state$response,
    delta1 = state$delta1,
    delta0 = state$delta0,
    par = state$par,
    alpha_patterns = state$alpha.patterns,
    zeta_patterns = state$zeta.patterns,
    patterns = state$patterns,
    step = theta_grid_length,
    lower = theta_lower,
    upper = theta_upper
  )
}

fcdcm_istem_update_parameters <- function(state) {
  sampled <- cpp_fcdcm_sample_alpha(
    theta = state$theta,
    response = state$response,
    delta1 = state$delta1,
    delta0 = state$delta0,
    alpha_patterns = state$alpha.patterns,
    zeta_patterns = state$zeta.patterns,
    class_prob = state$class.prob
  )
  state$alpha <- sampled$alpha
  state$zeta <- sampled$zeta
  state$class.index <- sampled$class.index

  delta <- fcdcm_istem_update_delta(
    theta = state$theta[, 1L],
    alpha = state$alpha,
    weight = istem_state_weight(state),
    delta1 = state$delta1,
    delta0 = state$delta0,
    prior = state$prior,
    optim.maxit = state$optim.maxit
  )
  state$delta1 <- delta$delta1
  state$delta0 <- delta$delta0

  state$par <- fcdcm_istem_update_par(
    response = state$response,
    weight = istem_state_weight(state),
    zeta = state$zeta,
    par = state$par,
    patterns = state$patterns,
    prior = state$prior
  )
  state$class.prob <- model.FCDCM(state$zeta.patterns, state$par,
                                  state$patterns)
  state
}

fcdcm_istem_update_delta <- function(theta, alpha, weight, delta1, delta0,
                                     prior, optim.maxit) {
  theta <- as.numeric(theta)
  alpha <- as.matrix(alpha)
  storage.mode(alpha) <- "integer"
  D <- ncol(alpha)

  for (d in seq_len(D)) {
    y <- alpha[, d]
    x0 <- c(log(pmax(delta1[d], 1e-8)), delta0[d])
    obj <- function(x) {
      if (any(!is.finite(x)) || any(abs(x) > 30)) return(Inf)
      slope <- exp(x[1L])
      location <- x[2L]
      eta <- slope * (theta - location)
      ll <- sum(weight * (y * stats::plogis(eta, log.p = TRUE) +
                            (1L - y) * stats::plogis(-eta, log.p = TRUE)))
      lp <- if (isTRUE(prior$use.prior)) {
        stats::dlnorm(slope, meanlog = prior$delta1.mu,
                      sdlog = prior$delta1.sigma, log = TRUE) +
          stats::dnorm(location, mean = prior$delta0.mu,
                       sd = prior$delta0.sigma, log = TRUE)
      } else {
        0
      }
      if (!is.finite(lp)) return(Inf)
      -ll - lp
    }
    opt <- tryCatch(
      stats::optim(
        par = x0,
        fn = obj,
        method = "BFGS",
        control = list(maxit = optim.maxit)
      ),
      error = function(e) NULL
    )
    if (!is.null(opt) && is.finite(opt$value)) {
      delta1[d] <- exp(opt$par[1L])
      delta0[d] <- opt$par[2L]
    }
  }

  list(delta1 = delta1, delta0 = delta0)
}

fcdcm_istem_update_par <- function(response, weight, zeta, par, patterns,
                                   prior) {
  response <- as.matrix(response)
  zeta <- as.matrix(zeta)
  storage.mode(response) <- "integer"
  storage.mode(zeta) <- "integer"
  patterns <- as.matrix(patterns)
  storage.mode(patterns) <- "integer"
  par <- fcdcm_normalize_par(par, B = nrow(patterns), arg = "par")

  for (b in seq_len(nrow(patterns))) {
    za <- zeta[, patterns[b, 1L]]
    zb <- zeta[, patterns[b, 2L]]
    y <- response[, b]
    low <- za < zb
    high <- za > zb
    par[b, "eta0"] <- fcdcm_istem_fit_eta(
      y = y[low],
      weight = weight[low],
      shift = 0,
      init = par[b, "eta0"],
      prior.mu = prior$eta0.mu,
      prior.sigma = prior$eta0.sigma,
      use.prior = prior$use.prior
    )
    par[b, "etaAB"] <- fcdcm_istem_fit_eta(
      y = y[high],
      weight = weight[high],
      shift = 0.5,
      init = par[b, "etaAB"],
      prior.mu = prior$etaAB.mu,
      prior.sigma = prior$etaAB.sigma,
      use.prior = prior$use.prior
    )
  }

  par
}

fcdcm_istem_fit_eta <- function(y, weight, shift, init, prior.mu,
                                prior.sigma, use.prior, upper = 0.5) {
  y <- as.integer(y)
  weight <- as.numeric(weight)
  eps <- 1e-6
  if (length(y) == 0L) {
    if (isTRUE(use.prior)) {
      return(max(eps, min(upper - eps,
               stats::rnorm(1L, mean = prior.mu, sd = prior.sigma))))
    }
    return(pmin(pmax(init, eps), upper - eps))
  }

  obj <- function(eta) {
    p <- pmin(pmax(shift + eta, eps), 1 - eps)
    ll <- sum(weight * (y * log(p) + (1L - y) * log1p(-p)))
    lp <- if (isTRUE(use.prior)) {
      stats::dnorm(eta, mean = prior.mu, sd = prior.sigma, log = TRUE)
    } else {
      0
    }
    if (!is.finite(lp)) return(Inf)
    -ll - lp
  }

  opt <- tryCatch(
    stats::optimize(obj, interval = c(eps, upper - eps)),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$objective)) {
    return(pmin(pmax(init, eps), upper - eps))
  }
  opt$minimum
}

fcdcm_istem_param_vec <- function(state) {
  c(
    as.numeric(state$delta1),
    as.numeric(state$delta0),
    as.numeric(state$par)
  )
}

fcdcm_istem_loglik_trace <- function(state, L, theta.lower, theta.upper) {
  theta.grid <- matrix(seq(theta.lower, theta.upper, length.out = L),
                       ncol = 1L)
  pi <- stats::dnorm(theta.grid[, 1L])
  pi <- pi / sum(pi)
  support <- fcdcm_latent_support(
    theta = theta.grid,
    pi = pi,
    delta1 = state$delta1,
    delta0 = state$delta0,
    par = state$par,
    alpha.patterns = state$alpha.patterns,
    zeta.patterns = state$zeta.patterns,
    patterns = state$patterns
  )
  istem_loglik_binary_sum(
    support$prob, state$response, state$response.count, support$pi
  )
}

fcdcm_istem_param_unpack <- function(pv, pv_se, pv_rhat = NULL, D, B,
                                     dim.names = NULL,
                                     block.names = NULL,
                                     par.names = c("eta0", "etaAB")) {
  idx <- 1L
  delta1 <- pv[idx:(idx + D - 1L)]
  delta1.se <- pv_se[idx:(idx + D - 1L)]
  delta1.rhat <- if (!is.null(pv_rhat)) pv_rhat[idx:(idx + D - 1L)] else rep(NA_real_, D)
  idx <- idx + D
  delta0 <- pv[idx:(idx + D - 1L)]
  delta0.se <- pv_se[idx:(idx + D - 1L)]
  delta0.rhat <- if (!is.null(pv_rhat)) pv_rhat[idx:(idx + D - 1L)] else rep(NA_real_, D)
  idx <- idx + D
  K <- length(par.names)
  par <- matrix(pv[idx:(idx + B * K - 1L)], nrow = B, ncol = K,
                dimnames = list(block.names, par.names))
  par.se <- matrix(pv_se[idx:(idx + B * K - 1L)], nrow = B, ncol = K,
                   dimnames = list(block.names, par.names))
  par.rhat <- if (!is.null(pv_rhat)) {
    matrix(pv_rhat[idx:(idx + B * K - 1L)], nrow = B, ncol = K,
           dimnames = list(block.names, par.names))
  } else {
    matrix(NA_real_, B, K, dimnames = list(block.names, par.names))
  }
  names(delta1) <- names(delta0) <- dim.names
  names(delta1.se) <- names(delta0.se) <- dim.names
  names(delta1.rhat) <- names(delta0.rhat) <- dim.names

  delta <- cbind(delta1 = delta1, delta0 = delta0)
  delta.se <- cbind(delta1 = delta1.se, delta0 = delta0.se)
  delta.Rhat <- cbind(delta1 = delta1.rhat, delta0 = delta0.rhat)
  list(
    delta = delta,
    delta.se = delta.se,
    delta.Rhat = delta.Rhat,
    par = par,
    par.se = par.se,
    par.Rhat = par.rhat
  )
}

fcdcm_posterior_alpha <- function(theta, response, delta1, delta0, par,
                                  alpha.patterns, zeta.patterns, patterns) {
  class.post <- fcdcm_class_posterior(
    theta = theta,
    response = response,
    delta1 = delta1,
    delta0 = delta0,
    par = par,
    alpha.patterns = alpha.patterns,
    zeta.patterns = zeta.patterns,
    patterns = patterns
  )
  alpha.prob <- class.post %*% as.matrix(alpha.patterns)
  colnames(alpha.prob) <- colnames(alpha.patterns)
  rownames(alpha.prob) <- rownames(response)
  alpha.prob
}

fcdcm_alpha_map <- function(class.post, alpha.patterns) {
  alpha.patterns <- as.matrix(alpha.patterns)
  alpha.est <- alpha.patterns[max.col(as.matrix(class.post),
                                      ties.method = "first"), , drop = FALSE]
  storage.mode(alpha.est) <- "integer"
  colnames(alpha.est) <- colnames(alpha.patterns)
  rownames(alpha.est) <- rownames(class.post)
  alpha.est
}

fcdcm_marginal_posterior <- function(response, delta1, delta0, par,
                                     alpha.patterns, zeta.patterns, patterns,
                                     L = 61, theta.lower = -6,
                                     theta.upper = 6, eps = 1e-12) {
  response <- as.matrix(response)
  storage.mode(response) <- "numeric"
  response.group <- istem_response_groups(response)
  y <- response.group$response
  count <- response.group$count
  G <- nrow(y)
  C <- nrow(alpha.patterns)

  theta.grid <- seq(theta.lower, theta.upper, length.out = L)
  pi.theta <- stats::dnorm(theta.grid)
  pi.theta <- pi.theta / sum(pi.theta)

  class.prob <- model.FCDCM(zeta.patterns, par, patterns)
  class.prob <- pmin(pmax(class.prob, eps), 1 - eps)
  log.p1 <- log(class.prob)
  log.p0 <- log1p(-class.prob)
  log.alpha <- fcdcm_alpha_profile_logprob(
    theta.grid, delta1, delta0, alpha.patterns, eps = eps
  )
  log.prior <- sweep(log.alpha, 1L, log(pi.theta), "+")

  class.post.group <- matrix(NA_real_, G, C)
  theta.group <- theta2.group <- numeric(G)

  for (g in seq_len(G)) {
    response.lp <- as.numeric(y[g, ] %*% t(log.p1) +
                                (1 - y[g, ]) %*% t(log.p0))
    log.post <- sweep(log.prior, 2L, response.lp, "+")
    max.log <- max(log.post)
    post <- exp(log.post - max.log)
    post <- post / sum(post)
    class.post.group[g, ] <- colSums(post)
    theta.weight <- rowSums(post)
    theta.group[g] <- sum(theta.grid * theta.weight)
    theta2.group[g] <- sum(theta.grid^2 * theta.weight)
  }

  class.post <- class.post.group[response.group$group, , drop = FALSE]
  colnames(class.post) <- pattern_key(alpha.patterns)
  rownames(class.post) <- rownames(response)

  theta <- matrix(theta.group[response.group$group], ncol = 1L)
  theta.se <- sqrt(pmax(theta2.group - theta.group^2, 0))
  theta.se <- matrix(theta.se[response.group$group], ncol = 1L)

  list(class.post = class.post, theta = theta, theta.se = theta.se,
       response.group = response.group, count = count)
}

fcdcm_class_posterior <- function(theta, response, delta1, delta0, par,
                                  alpha.patterns, zeta.patterns, patterns,
                                  eps = 1e-12) {
  class.prob <- model.FCDCM(zeta.patterns, par, patterns)
  fcdcm_class_posterior_class_prob(
    theta = theta,
    response = response,
    delta1 = delta1,
    delta0 = delta0,
    alpha.patterns = alpha.patterns,
    class.prob = class.prob,
    eps = eps
  )
}

fcdcm_class_posterior_class_prob <- function(theta, response, delta1, delta0,
                                             alpha.patterns, class.prob,
                                             eps = 1e-12) {
  theta <- as.numeric(theta)
  response <- as.matrix(response)
  storage.mode(response) <- "numeric"
  alpha.patterns <- as.matrix(alpha.patterns)
  class.prob <- as.matrix(class.prob)
  class.prob <- pmin(pmax(class.prob, eps), 1 - eps)

  log.alpha <- fcdcm_alpha_profile_logprob(
    theta, delta1, delta0, alpha.patterns, eps = eps
  )
  log.p1 <- log(class.prob)
  log.p0 <- log1p(-class.prob)
  log.response <- response %*% t(log.p1) + (1 - response) %*% t(log.p0)
  log.post <- log.alpha + log.response
  max.log <- apply(log.post, 1L, max)
  post <- exp(sweep(log.post, 1L, max.log, "-"))
  post <- post / rowSums(post)

  colnames(post) <- pattern_key(alpha.patterns)
  rownames(post) <- rownames(response)
  post
}

fcdcm_stan_class_posterior <- function(theta.draws, delta1.draws, delta0.draws,
                                       eta0.draws, etaAB.draws,
                                       response,
                                       alpha.patterns, zeta.patterns,
                                       patterns, eps = 1e-12) {
  theta.draws <- as.matrix(theta.draws)
  delta1.draws <- as.matrix(delta1.draws)
  delta0.draws <- as.matrix(delta0.draws)
  eta0.draws <- as.matrix(eta0.draws)
  etaAB.draws <- as.matrix(etaAB.draws)

  S <- nrow(theta.draws)
  response <- as.matrix(response)
  N <- nrow(response)
  C <- nrow(alpha.patterns)
  class.post <- matrix(0, N, C)

  for (s in seq_len(S)) {
    par.s <- cbind(eta0 = eta0.draws[s, ], etaAB = etaAB.draws[s, ])
    class.post <- class.post + fcdcm_class_posterior(
      theta = theta.draws[s, ],
      response = response,
      delta1 = delta1.draws[s, ],
      delta0 = delta0.draws[s, ],
      par = par.s,
      alpha.patterns = alpha.patterns,
      zeta.patterns = zeta.patterns,
      patterns = patterns,
      eps = eps
    )
  }

  class.post <- class.post / S
  colnames(class.post) <- pattern_key(alpha.patterns)
  rownames(class.post) <- rownames(response)
  class.post
}
