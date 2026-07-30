#' Fit FCDCM Model via Stan
#'
#' @noRd
fit.FCDCM.stan <- function(data, Q.matrix, block.items = NULL,
                           dcm.type = "DINA",
                           control.model = NULL,
                           control.method = NULL,
                           .call = NULL) {
  call <- if (is.null(.call)) match.call() else .call
  control <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  stan.method <- fc_stan_method_control(control.method)
  prior <- list(
    delta1.mu = istem_scalar(get_ctrl("delta1.mu", 0.0, control), "delta1.mu"),
    delta1.sigma = istem_scalar(get_ctrl("delta1.sigma", 0.50, control), "delta1.sigma", positive = TRUE),
    delta0.mu = istem_scalar(get_ctrl("delta0.mu", 0.00, control), "delta0.mu"),
    delta0.sigma = istem_scalar(get_ctrl("delta0.sigma", 1.00, control), "delta0.sigma", positive = TRUE),
    eta0.mu = istem_scalar(
      get_ctrl("eta0.mu", 0.10, control),
      "eta0.mu", positive = TRUE),
    eta0.sigma = istem_scalar(
      get_ctrl("eta0.sigma", 0.10, control),
      "eta0.sigma", positive = TRUE),
    etaAB.mu = istem_scalar(
      get_ctrl("etaAB.mu", 0.30, control),
      "etaAB.mu", positive = TRUE),
    etaAB.sigma = istem_scalar(
      get_ctrl("etaAB.sigma", 0.10, control),
      "etaAB.sigma", positive = TRUE),
    eta_equal0.mu = istem_scalar(
      get_ctrl("eta_equal0.mu", 0.50, control),
      "eta_equal0.mu", positive = TRUE),
    eta_equal0.sigma = istem_scalar(
      get_ctrl("eta_equal0.sigma", 0.25, control),
      "eta_equal0.sigma", positive = TRUE),
    eta_equal1.mu = istem_scalar(
      get_ctrl("eta_equal1.mu", 0.50, control),
      "eta_equal1.mu", positive = TRUE),
    eta_equal1.sigma = istem_scalar(
      get_ctrl("eta_equal1.sigma", 0.25, control),
      "eta_equal1.sigma", positive = TRUE)
  )

  sm <- stanmodels[["FCDCM"]]
  if (is.null(sm)) {
    stop("Stan model FCDCM not found !", call. = FALSE)
  }

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
  mcmc_env <- setup_mcmc_env(common.method$cores)
  on.exit(restore_mcmc_env(mcmc_env), add = TRUE)

  stan.control <- build_stan_control(stan.method$algorithm, control.method)

  stan.data <- list(
    N = prep$N,
    B = prep$B,
    I = prep$I.states,
    D = prep$D,
    C = nrow(prep$alpha.patterns),
    y = as.integer(prep$response),
    alpha_patterns = as.integer(t(prep$alpha.patterns)),
    zeta_patterns = as.integer(t(prep$zeta.patterns)),
    patterns = as.integer(t(prep$patterns)),
    delta1_mu = prior$delta1.mu,
    delta1_sigma = prior$delta1.sigma,
    delta0_mu = prior$delta0.mu,
    delta0_sigma = prior$delta0.sigma,
    eta0_mu = prior$eta0.mu,
    eta0_sigma = prior$eta0.sigma,
    etaAB_mu = prior$etaAB.mu,
    etaAB_sigma = prior$etaAB.sigma,
    eta_equal0_mu = prior$eta_equal0.mu,
    eta_equal0_sigma = prior$eta_equal0.sigma,
    eta_equal1_mu = prior$eta_equal1.mu,
    eta_equal1_sigma = prior$eta_equal1.sigma,
    free_equal = 0L
  )

  set.seed(common.method$seed)
  stan.obj <- suppressWarnings(
    rstan::sampling(
      object = sm,
      data = stan.data,
      chains = stan.method$chains,
      iter = stan.method$iter,
      warmup = stan.method$warmup,
      thin = stan.method$thin,
      init = stan.method$init,
      algorithm = stan.method$algorithm,
      cores = mcmc_env$cores,
      control = stan.control,
      verbose = common.method$vis,
      show_messages = common.method$vis,
      seed = common.method$seed
    )
  )

  MCMC.obj <- rstan::extract(stan.obj, permuted = TRUE)
  stan.sum <- rstan::summary(stan.obj)$summary

  theta.draws <- as.matrix(MCMC.obj$theta)
  theta <- matrix(colMeans(theta.draws), ncol = 1L)
  theta.se <- matrix(apply(theta.draws, 2L, stats::sd), ncol = 1L)
  theta.Rhat <- matrix(extract_rhat_vector(stan.sum, "theta", prep$N),
                       ncol = 1L)
  colnames(theta) <- colnames(theta.se) <- colnames(theta.Rhat) <- "theta"
  rownames(theta) <- rownames(theta.se) <- rownames(theta.Rhat) <-
    rownames(prep$response)

  delta1.draws <- as.matrix(MCMC.obj$delta1)
  delta0.draws <- as.matrix(MCMC.obj$delta0)
  delta <- cbind(
    delta1 = colMeans(delta1.draws),
    delta0 = colMeans(delta0.draws)
  )
  delta.se <- cbind(
    delta1 = apply(delta1.draws, 2L, stats::sd),
    delta0 = apply(delta0.draws, 2L, stats::sd)
  )
  delta.Rhat <- cbind(
    delta1 = extract_rhat_vector(stan.sum, "delta1", prep$D),
    delta0 = extract_rhat_vector(stan.sum, "delta0", prep$D)
  )
  rownames(delta) <- rownames(delta.se) <- rownames(delta.Rhat) <-
    prep$dim.names

  eta0.draws <- as.matrix(MCMC.obj$eta0)
  etaAB.draws <- as.matrix(MCMC.obj$etaAB)
  par <- cbind(eta0 = colMeans(eta0.draws),
               etaAB = colMeans(etaAB.draws))
  par.se <- cbind(eta0 = apply(eta0.draws, 2L, stats::sd),
                  etaAB = apply(etaAB.draws, 2L, stats::sd))
  par.Rhat <- cbind(eta0 = extract_rhat_vector(stan.sum, "eta0", prep$B),
                    etaAB = extract_rhat_vector(stan.sum, "etaAB", prep$B))
  rownames(par) <- rownames(par.se) <- rownames(par.Rhat) <-
    prep$block.names
  par.free <- matrix(TRUE, prep$B, ncol(par),
                     dimnames = list(prep$block.names,
                                     colnames(par)))

  class.post <- fcdcm_stan_class_posterior(
    theta.draws = theta.draws,
    delta1.draws = delta1.draws,
    delta0.draws = delta0.draws,
    eta0.draws = eta0.draws,
    etaAB.draws = etaAB.draws,
    response = prep$response,
    alpha.patterns = prep$alpha.patterns,
    zeta.patterns = prep$zeta.patterns,
    patterns = prep$patterns
  )
  colnames(class.post) <- pattern_key(prep$alpha.patterns)
  rownames(class.post) <- rownames(prep$response)
  alpha.prob <- class.post %*% prep$alpha.patterns
  colnames(alpha.prob) <- prep$dim.names
  rownames(alpha.prob) <- rownames(prep$response)
  alpha.est <- fcdcm_alpha_map(class.post, prep$alpha.patterns)
  storage.mode(alpha.est) <- "integer"
  alpha.se <- matrix(NA_real_, prep$N, prep$D,
                     dimnames = dimnames(alpha.est))
  alpha.Rhat <- alpha.se

  npar <- ncol(par) * prep$B + 2L * prep$D
  results <- list(
    npar = npar,
    method = "stan",
    theta = list(est = theta, se = theta.se, Rhat = theta.Rhat),
    alpha = list(est = alpha.est, se = alpha.se, Rhat = alpha.Rhat,
                 prob = alpha.prob),
    class.post = class.post,
    delta = list(est = delta, se = delta.se, Rhat = delta.Rhat),
    par = list(est = par, se = par.se, Rhat = par.Rhat, free = par.free),
    Corr = list(
      est = matrix(1, 1L, 1L, dimnames = list("theta", "theta")),
      se = matrix(0, 1L, 1L, dimnames = list("theta", "theta")),
      Rhat = matrix(NA_real_, 1L, 1L, dimnames = list("theta", "theta"))
    ),
    stan.obj = stan.obj,
    MCMC.obj = MCMC.obj,
    log_lik = MCMC.obj$log_lik,
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
      method = "stan",
      chains = stan.method$chains,
      iter = stan.method$iter,
      warmup = stan.method$warmup,
      thin = stan.method$thin,
      init = stan.method$init,
      algorithm = stan.method$algorithm,
      cores = common.method$cores,
      vis = common.method$vis,
      seed = common.method$seed,
      control.model = control,
      control.method = fit_effective_method_control(
        control.method, method = stan.method, common = common.method
      )
    )
  )

  class(results) <- "FCDCM"
  L <- get_ctrl("L", NULL, control)
  results$logLik <- logLik.FCDCM(results, L = L)
  results
}
