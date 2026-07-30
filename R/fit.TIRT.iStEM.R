################################# TIRT iStEM #################################

fit.TIRT.iStEM <- function(data, Q.matrix, block.items = NULL,
                           fc.type = "RANK",
                           control.model = NULL, control.method = NULL,
                           .call = NULL) {

  call <- if (is.null(.call)) match.call() else .call
  control.model <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  vis <- common.method$vis
  method <- istem_method_control(control.method)
  set.seed(common.method$seed)

  prep <- tirt_istem_prepare_data(data, Q.matrix, block.items, fc.type,
                                  control.model)
  D <- prep$D
  N <- prep$N.person
  G <- prep$G
  I.states <- prep$I.states
  method <- istem_apply_grid_control(method, control.model, D)
  istem_check_method_control(method)
  control.model <- istem_effective_grid_model_control(control.model, method)

  prior <- tirt_istem_prior_control(control.model)
  bounds <- list(
    lambda.lower = get_ctrl("lambda.lower", prior$lambda.alpha, control.method),
    lambda.upper = get_ctrl("lambda.upper", prior$lambda.beta, control.method),
    psi.lower = get_ctrl("psi.lower", 0.01, control.method),
    psi.upper = get_ctrl("psi.upper", 5, control.method),
    gamma.lower = get_ctrl("gamma.lower", -6, control.method),
    gamma.upper = get_ctrl("gamma.upper", 6, control.method)
  )
  tirt_istem_check_bounds(bounds)

  theta_mu <- array(get_ctrl("theta.mu", rep(0, D), control.model), dim = D)
  init <- tirt_istem_initial_par(prep, prior, bounds)
  include_corr <- D > 1L && !method$fix.corr

  state <- c(
    prep,
    init,
    list(
      theta = istem_random_theta(G, D, theta_mu, method$theta.lower,
                                 method$theta.upper),
      Corr = istem_random_corr(D, include_corr),
      theta_mu = theta_mu,
      optim.maxit = method$optim.maxit,
      bounds = bounds,
      prior = prior
    )
  )

  run <- istem_run(
    state = state,
    method = method,
    N = N,
    vis = vis,
    label = "TIRT",
    update_parameters = tirt_istem_update_parameters,
    sample_theta = tirt_istem_sample_theta,
    param_vec = function(state) {
      tirt_istem_param_vec(state$lambda, state$psi, state$gamma,
                           state$Idx.lambda.est, state$Idx.psi.est,
                           state$Idx.gamma.est, state$Corr, include_corr)
    },
    logLik_fun = function(state) {
      tirt_istem_loglik_trace(
        state, method$L, method$theta.lower, method$theta.upper
      )
    }
  )

  final <- tirt_istem_param_unpack(
    pv = run$chain.mean,
    pv_se = run$chain.sd,
    pv_rhat = NULL,
    state = run$state,
    include_corr = include_corr
  )

  theta.est <- istem_expand_group_matrix(
    run$theta.est, prep$response.group, rownames(prep$response))
  theta.se <- istem_expand_group_matrix(
    run$theta.se, prep$response.group, rownames(prep$response))
  theta_est <- istem_named_theta(theta.est, theta.se, prep$response, D)

  lambda_sign <- apply(prep$Q.matrix, 1, function(x) {
    nz <- which(x != 0)
    if (length(nz) == 0L) return(1)
    sign(x[nz[1L]])
  })
  lambda_signed <- final$lambda * lambda_sign
  lambda.free <- rep(FALSE, I.states)
  lambda.free[prep$Idx.lambda.est] <- TRUE
  psi.free <- rep(FALSE, I.states)
  psi.free[prep$Idx.psi.est] <- TRUE
  gamma.free <- rep(FALSE, length(final$gamma))
  gamma.free[prep$Idx.gamma.est] <- TRUE

  psi2 <- final$psi^2
  psi2.se <- 2 * final$psi * final$psi.se
  psi2.rhat <- final$psi.rhat
  final$lambda.se[!lambda.free] <- NA_real_
  final$lambda.rhat[!lambda.free] <- NA_real_
  psi2.se[!psi.free] <- NA_real_
  psi2.rhat[!psi.free] <- NA_real_
  final$gamma.se[!gamma.free] <- NA_real_
  final$gamma.rhat[!gamma.free] <- NA_real_

  par_est <- list(
    est = cbind(lambda_signed * abs(prep$Q.matrix), psi2),
    se = cbind(final$lambda.se * abs(prep$Q.matrix), psi2.se),
    Rhat = cbind(final$lambda.rhat * abs(prep$Q.matrix), psi2.rhat),
    free = tirt_par_free_mask(prep$Q.matrix, lambda.free, psi.free)
  )
  colnames(par_est$est) <- colnames(par_est$se) <- colnames(par_est$Rhat) <-
    colnames(par_est$free) <-
    c(paste0("Dim.", seq_len(D)), "psi2")
  rownames(par_est$est) <- rownames(par_est$se) <- rownames(par_est$Rhat) <-
    rownames(par_est$free) <- paste0("item", seq_len(I.states))
  par_est$se[!par_est$free] <- NA_real_
  par_est$Rhat[!par_est$free] <- NA_real_

  gamma.matrix <- tirt_istem_gamma_matrix(final$gamma, final$gamma.se,
                                          prep$pairs.matrix, I.states)

  Corr.Rhat <- final$Corr.rhat
  dimnames(Corr.Rhat) <- list(paste0("Dim.", seq_len(D)), paste0("Dim.", seq_len(D)))

  results <- list(
    npar = length(run$chain.mean),
    method = "iStEM",
    theta = theta_est,
    par = par_est,
    gamma = list(est = final$gamma, se = final$gamma.se,
                 Rhat = final$gamma.rhat,
                 free = gamma.free),
    Corr = list(est = final$Corr, se = final$Corr.se, Rhat = Corr.Rhat),
    gamma.matrix = list(est = gamma.matrix$est, se = gamma.matrix$se,
                        Rhat = gamma.matrix$Rhat),
    Q.matrix = prep$Q.matrix,
    block.items = prep$block.items,
    pairs.matrix = prep$pairs.matrix,
    pairs.value = prep$pairs.value,
    response = prep$response,
    fc.type = prep$fc.type,
    stan.obj = NULL,
    MCMC.obj = NULL,
    call = call,
    arguments = list(
      data = data,
      Q.matrix = prep$Q.matrix,
      block.items = prep$block.items,
      fc.type = prep$fc.type,
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
  class(results) <- "TIRT"

  results$logLik <- logLik(results, L = method$L,
                           theta.low = method$theta.lower,
                           theta.up = method$theta.upper)
  results
}

tirt_istem_prior_control <- function(control.model) {
  lambda_alpha <- istem_scalar(get_ctrl("lambda.alpha", 0.40, control.model),
                               "lambda.alpha", positive = TRUE)
  lambda_beta <- istem_scalar(get_ctrl("lambda.beta", 0.90, control.model),
                              "lambda.beta", positive = TRUE)
  if (lambda_alpha >= lambda_beta) {
    stop("'lambda.alpha' must be smaller than 'lambda.beta'.",
         call. = FALSE)
  }
  list(
    use.prior = isTRUE(get_ctrl("use.prior", TRUE, control.model)),
    lambda.alpha = lambda_alpha,
    lambda.beta = lambda_beta,
    psi.mu = istem_scalar(get_ctrl("psi.mu", 0.80, control.model), "psi.mu",
                          positive = TRUE),
    psi.sigma = istem_scalar(get_ctrl("psi.sigma", 0.30, control.model),
                             "psi.sigma", positive = TRUE),
    gamma.mu = istem_scalar(get_ctrl("gamma.mu", 0, control.model),
                            "gamma.mu"),
    gamma.sigma = istem_scalar(get_ctrl("gamma.sigma", 0.8, control.model),
                               "gamma.sigma", positive = TRUE)
  )
}

tirt_istem_check_bounds <- function(bounds) {
  for (nm in names(bounds)) {
    bounds[[nm]] <- istem_scalar(bounds[[nm]], nm, arg = "control.method")
  }
  if (bounds$lambda.lower <= 0 ||
      bounds$lambda.lower >= bounds$lambda.upper) {
    stop("'lambda.lower' must be positive and smaller than 'lambda.upper'.",
         call. = FALSE)
  }
  if (bounds$psi.lower <= 0 || bounds$psi.lower >= bounds$psi.upper) {
    stop("'psi.lower' must be positive and smaller than 'psi.upper'.",
         call. = FALSE)
  }
  if (bounds$gamma.lower >= bounds$gamma.upper) {
    stop("'gamma.lower' must be smaller than 'gamma.upper'.", call. = FALSE)
  }
  invisible(NULL)
}

tirt_istem_prepare_data <- function(data, Q.matrix, block.items, fc.type,
                                    control.model) {
  fc.type <- toupper(fc.type)
  N.block.guess <- if (is.list(block.items)) length(block.items) else ncol(data)
  if (length(fc.type) != N.block.guess) {
    fc.type <- rep(fc.type[1L], N.block.guess)
  }
  if (is.null(block.items)) {
    block.items <- get.block.items.from.data(data)
  }
  result <- get.response.from.data.TIRT(data, block.items, fc.type = fc.type)
  response <- result$response
  pairs.value <- result$pairs.value

  if (any(is.na(response))) {
    first_missing <- which(is.na(response), arr.ind = TRUE)[1L, ]
    stop(
      sprintf(
        "TIRT responses must not contain missing values (row %d, column %d).",
        first_missing[1L], first_missing[2L]
      ),
      call. = FALSE
    )
  }
  if (!all(response %in% c(0, 1))) {
    stop("TIRT responses must contain only 0 and 1.",
         call. = FALSE)
  }

  if (!is.matrix(Q.matrix)) {
    stop("'Q.matrix' must be a matrix.", call. = FALSE)
  }
  if (!all(Q.matrix %in% c(-1, 0, 1))) {
    stop("'Q.matrix' must contain only -1, 0, or 1.", call. = FALSE)
  }

  I.states <- nrow(Q.matrix)
  D <- ncol(Q.matrix)
  if (I.states < 1L) {
    stop("'Q.matrix' must have at least one row (item).", call. = FALSE)
  }
  if (D < 2L) {
    stop("At least two trait dimensions must be specified in Q.matrix (TIRT requires D >= 2).",
         call. = FALSE)
  }
  if (any(rowSums(abs(Q.matrix)) < 1L)) {
    stop("Each item must measure at least one trait (row sums of abs('Q.matrix') must be >= 1).",
         call. = FALSE)
  }

  if (!is.list(block.items) || length(block.items) < 1L) {
    stop("'block.items' must be a non-empty list.", call. = FALSE)
  }
  if (!all(vapply(block.items, function(x) length(x) >= 2L, logical(1L)))) {
    stop("Each block in 'block.items' must contain at least 2 items.",
         call. = FALSE)
  }
  all_items_in_blocks <- unlist(block.items, use.names = FALSE)
  if (!all(all_items_in_blocks %in% seq_len(I.states))) {
    stop("Item indices in 'block.items' must be within 1:I.states.",
         call. = FALSE)
  }
  if (anyDuplicated(all_items_in_blocks)) {
    stop("Item indices in 'block.items' must be unique across all blocks.",
         call. = FALSE)
  }
  if (length(all_items_in_blocks) != I.states) {
    stop("The total number of items in 'block.items' must exactly match nrow(Q.matrix).",
         call. = FALSE)
  }

  I.block <- length(block.items[[1L]])
  if (!all(vapply(block.items, length, integer(1L)) == I.block)) {
    stop("All blocks in 'block.items' must have the same number of items.",
         call. = FALSE)
  }
  N.block <- length(block.items)
  if (length(fc.type) != N.block) {
    fc.type <- rep(fc.type[1L], N.block)
  }

  I.block.vec <- vapply(block.items, length, integer(1L))
  block.pair.counts <- integer(N.block)
  for (b in seq_len(N.block)) {
    k <- I.block.vec[b]
    block.pair.counts[b] <- switch(fc.type[b],
      RANK = choose(k, 2),
      MOLE = 2L * k - 3L,
      PICK = k - 1L,
      stop("'fc.type' must be RANK, MOLE, or PICK.", call. = FALSE)
    )
  }
  N.pairs.theoretical <- sum(block.pair.counts)
  if (ncol(response) != N.pairs.theoretical) {
    stop("Column dimension of response matrix does not match block.items/fc.type.",
         call. = FALSE)
  }

  response.raw <- response
  pairs.value.raw <- pairs.value
  if (!is.null(pairs.value) && !all(fc.type == "RANK")) {
    response.group <- istem_response_pair_groups(response, pairs.value)
    pairs.value <- response.group$pairs.value
  } else {
    response.group <- istem_response_groups(response)
    if (!is.null(pairs.value)) {
      pairs.value <- pairs.value[response.group$first]
    }
  }
  response <- response.group$response

  N.pairs.all <- sum(choose(I.block.vec, 2))
  pairs.matrix <- matrix(0L, N.pairs.all, 2L)
  idx <- 0L
  for (b in seq_len(N.block)) {
    pairs.cur <- t(combn(block.items[[b]], 2))
    pairs.matrix[idx + seq_len(nrow(pairs.cur)), ] <- pairs.cur
    idx <- idx + nrow(pairs.cur)
  }
  pair.keys <- paste0(pairs.matrix[, 1L], "_", pairs.matrix[, 2L])

  N.person <- nrow(response.raw)
  G <- response.group$G
  N.obs <- G * ncol(response)
  Y <- integer(N.obs)
  Idx.person <- integer(N.obs)
  Idx.item.i <- integer(N.obs)
  Idx.item.k <- integer(N.obs)
  Idx.pair <- integer(N.obs)

  idx <- 0L
  col_offset <- 0L
  for (b in seq_len(N.block)) {
    n_pairs.b <- block.pair.counts[b]
    if (fc.type[b] == "RANK") {
      pairs.b <- t(combn(block.items[[b]], 2))
      for (p in seq_len(G)) {
        for (pp in seq_len(n_pairs.b)) {
          idx <- idx + 1L
          Y[idx] <- response[p, col_offset + pp]
          Idx.person[idx] <- p
          Idx.item.i[idx] <- pairs.b[pp, 1L]
          Idx.item.k[idx] <- pairs.b[pp, 2L]
          Idx.pair[idx] <- match(paste0(pairs.b[pp, 1L], "_",
                                        pairs.b[pp, 2L]), pair.keys)
        }
      }
    } else {
      for (p in seq_len(G)) {
        pairs.cur <- pairs.value[[p]][[b]]
        if (!is.matrix(pairs.cur)) pairs.cur <- matrix(pairs.cur, ncol = 2)
        for (pp in seq_len(n_pairs.b)) {
          idx <- idx + 1L
          Y[idx] <- response[p, col_offset + pp]
          Idx.person[idx] <- p
          Idx.item.i[idx] <- pairs.cur[pp, 1L]
          Idx.item.k[idx] <- pairs.cur[pp, 2L]
          Idx.pair[idx] <- match(paste0(pairs.cur[pp, 1L], "_",
                                        pairs.cur[pp, 2L]), pair.keys)
        }
      }
    }
    col_offset <- col_offset + n_pairs.b
  }

  Idx.gamma.est <- sort(unique(Idx.pair))
  Idx.gamma.fix <- setdiff(seq_len(N.pairs.all), Idx.gamma.est)

  Idx.lambda.est.raw <- which(rowSums(abs(Q.matrix)) > 0)
  is.case.b <- D == 2L && I.block == 2L
  if (is.case.b) {
    first_pair_items <- block.items[[1L]]
    Idx.lambda.fix <- first_pair_items
    lambda.fix.val <- rep(0.80, length(first_pair_items))
    Idx.lambda.est <- setdiff(Idx.lambda.est.raw, Idx.lambda.fix)
  } else {
    Idx.lambda.est <- Idx.lambda.est.raw
    Idx.lambda.fix <- integer(0)
    lambda.fix.val <- numeric(0)
  }

  last_item_iblock <- unname(vapply(block.items, function(x) x[length(x)],
                                    integer(1L)))
  is.case.a <- I.block == 2L && D > 2L
  if (is.case.a || is.case.b) {
    Idx.psi.fix <- seq_len(I.states)
    psi.fix.val <- rep(sqrt(0.5), I.states)
    Idx.psi.est <- integer(0)
  } else {
    Idx.psi.fix <- last_item_iblock
    psi.fix.val <- rep(1, length(last_item_iblock))
    Idx.psi.est <- setdiff(seq_len(I.states), Idx.psi.fix)
  }

  gamma.mu <- get_ctrl("gamma.mu", 0, control.model)
  gamma.fix.val <- rep(gamma.mu, length(Idx.gamma.fix))
  if (length(Idx.gamma.fix) > 0L) {
    warning(sprintf(
      "TIRT: %d of %d possible pairwise comparisons were never observed. The corresponding gamma parameters are fixed to the prior mean (%.2f).",
      length(Idx.gamma.fix), N.pairs.all, gamma.mu
    ))
  }

  obs.by.person <- split(seq_along(Y), Idx.person)

  list(
    response = response.raw,
    response.unique = response,
    response.count = response.group$count,
    response.group = response.group$group,
    G = G,
    Q.matrix = Q.matrix,
    block.items = block.items,
    pairs.value = pairs.value.raw,
    pairs.value.unique = pairs.value,
    pairs.matrix = pairs.matrix,
    fc.type = fc.type,
    N.person = N.person,
    I.states = I.states,
    D = D,
    Y = Y,
    Idx.person = Idx.person,
    Idx.item.i = Idx.item.i,
    Idx.item.k = Idx.item.k,
    Idx.pair = Idx.pair,
    obs.by.person = obs.by.person,
    Idx.lambda.est = Idx.lambda.est,
    Idx.lambda.fix = Idx.lambda.fix,
    lambda.fix.val = lambda.fix.val,
    Idx.psi.est = Idx.psi.est,
    Idx.psi.fix = Idx.psi.fix,
    psi.fix.val = psi.fix.val,
    Idx.gamma.est = Idx.gamma.est,
    Idx.gamma.fix = Idx.gamma.fix,
    gamma.fix.val = gamma.fix.val
  )
}

tirt_istem_initial_par <- function(prep, prior, bounds) {
  lambda <- stats::runif(prep$I.states, bounds$lambda.lower,
                         bounds$lambda.upper)
  if (length(prep$Idx.lambda.fix) > 0L) {
    lambda[prep$Idx.lambda.fix] <- prep$lambda.fix.val
  }
  psi <- istem_rnorm_bounded(prep$I.states, prior$psi.mu, prior$psi.sigma,
                             bounds$psi.lower, bounds$psi.upper)
  if (length(prep$Idx.psi.fix) > 0L) {
    psi[prep$Idx.psi.fix] <- prep$psi.fix.val
  }
  gamma <- istem_rnorm_bounded(nrow(prep$pairs.matrix), prior$gamma.mu,
                               prior$gamma.sigma, bounds$gamma.lower,
                               bounds$gamma.upper)
  if (length(prep$Idx.gamma.fix) > 0L) {
    gamma[prep$Idx.gamma.fix] <- prep$gamma.fix.val
  }
  list(lambda = lambda, psi = psi, gamma = gamma)
}

tirt_istem_sample_theta <- function(state, theta_grid_length,
                                    theta_lower, theta_upper) {
  chol <- istem_corr_chol(state$Corr)
  cpp_gibbs_tirt_theta(
    theta_ = state$theta,
    Q_matrix_ = state$Q.matrix,
    lambda_ = state$lambda,
    psi_ = state$psi,
    gamma_ = state$gamma,
    Y_ = state$Y,
    obs_by_person_ = state$obs.by.person,
    Idx_item_i_ = state$Idx.item.i,
    Idx_item_k_ = state$Idx.item.k,
    Idx_pair_ = state$Idx.pair,
    theta_mu_ = state$theta_mu,
    chol_corr_ = chol$chol,
    log_diag_sum = chol$log_diag_sum,
    step = theta_grid_length,
    lower = theta_lower,
    upper = theta_upper
  )
}

tirt_istem_update_parameters <- function(state) {
  # Pre-compute qtheta once (changes only when theta changes, i.e. between batches)
  qtheta <- state$theta %*% t(state$Q.matrix)

  items_to_opt <- union(state$Idx.lambda.est, state$Idx.psi.est)

  # ---- Build static observation caches (invariant across coord-ascent passes) ----
  # Per-item cache: for each item, the observation indices and structural info
  item_cache <- vector("list", length(items_to_opt))
  names(item_cache) <- as.character(items_to_opt)

  for (idx in seq_along(items_to_opt)) {
    i <- items_to_opt[idx]
    obs_idx <- which(state$Idx.item.i == i | state$Idx.item.k == i)
    if (length(obs_idx) == 0L) {
      item_cache[[idx]] <- NULL
      next
    }

    persons  <- state$Idx.person[obs_idx]
    ii       <- state$Idx.item.i[obs_idx]
    ik       <- state$Idx.item.k[obs_idx]
    ip       <- state$Idx.pair[obs_idx]
    is_first <- ii == i

    # qt_target and qt_other are invariant (theta is fixed during param update)
    qt_target <- numeric(length(obs_idx))
    qt_other  <- numeric(length(obs_idx))

    idx_f <- which(is_first)
    if (length(idx_f) > 0L) {
      qt_target[idx_f] <- qtheta[cbind(persons[idx_f], ii[idx_f])]
      qt_other[idx_f]  <- qtheta[cbind(persons[idx_f], ik[idx_f])]
    }
    idx_s <- which(!is_first)
    if (length(idx_s) > 0L) {
      qt_target[idx_s] <- qtheta[cbind(persons[idx_s], ik[idx_s])]
      qt_other[idx_s]  <- qtheta[cbind(persons[idx_s], ii[idx_s])]
    }

    item_cache[[idx]] <- list(
      i         = i,
      qt_target = qt_target,
      qt_other  = qt_other,
      other_item = ifelse(is_first, ik, ii),   # the OTHER item in each obs
      pair_idx   = ip,
      is_first   = as.integer(is_first),
      Y          = state$Y[obs_idx],
      weight     = state$response.count[persons],
      opt_lambda = i %in% state$Idx.lambda.est,
      opt_psi    = i %in% state$Idx.psi.est
    )
  }

  # Per-pair cache
  pair_cache <- vector("list", length(state$Idx.gamma.est))
  names(pair_cache) <- as.character(state$Idx.gamma.est)

  for (idx in seq_along(state$Idx.gamma.est)) {
    ip <- state$Idx.gamma.est[idx]
    obs_idx <- which(state$Idx.pair == ip)
    if (length(obs_idx) == 0L) {
      pair_cache[[idx]] <- NULL
      next
    }

    persons <- state$Idx.person[obs_idx]
    ii      <- state$Idx.item.i[obs_idx]
    ik      <- state$Idx.item.k[obs_idx]

    # Pre-compute the part of mu that is independent of gamma
    qt_i <- qtheta[cbind(persons, ii)]
    qt_k <- qtheta[cbind(persons, ik)]

    pair_cache[[idx]] <- list(
      ip      = ip,
      persons = persons,
      ii      = ii,
      ik      = ik,
      qt_i    = qt_i,
      qt_k    = qt_k,
      Y       = state$Y[obs_idx],
      weight  = state$response.count[persons]
    )
  }

  # ---- Coordinate ascent: iterate items -> pairs until convergence ----
  max_ca_passes <- 10L
  for (ca_pass in seq_len(max_ca_passes)) {
    max_change <- 0.0

    # -- Pass 1: update all item (lambda, psi) --
    for (cache in item_cache) {
      if (is.null(cache)) next
      i <- cache$i

      # Recompute parameter-dependent arrays (may have changed since last pass)
      lambda_other <- state$lambda[cache$other_item]
      psi_other    <- state$psi[cache$other_item]
      gamma_pair   <- state$gamma[cache$pair_idx]

      old_lambda <- state$lambda[i]
      old_psi    <- state$psi[i]

      if (cache$opt_lambda && cache$opt_psi) {
        x0    <- c(old_lambda, old_psi)
        lower <- c(state$bounds$lambda.lower, state$bounds$psi.lower)
        upper <- c(state$bounds$lambda.upper, state$bounds$psi.upper)

        obj <- function(x) {
          ll <- cpp_tirt_item_loglik_weighted(
            cache$qt_target, cache$qt_other,
            lambda_other, psi_other, gamma_pair,
            cache$is_first, cache$Y, cache$weight, x[1L], x[2L])
          lp <- if (isTRUE(state$prior$use.prior)) {
            stats::dnorm(x[2L], mean = state$prior$psi.mu,
                         sd = state$prior$psi.sigma, log = TRUE)
          } else 0
          if (!is.finite(lp)) return(Inf)
          -ll - lp
        }

        opt <- tryCatch(
          stats::optim(par = x0, fn = obj, method = "L-BFGS-B",
                       lower = lower, upper = upper,
                       control = list(maxit = state$optim.maxit)),
          error = function(e) NULL
        )
        if (!is.null(opt) && is.finite(opt$value)) {
          state$lambda[i] <- opt$par[1L]
          state$psi[i]    <- opt$par[2L]
        }

      } else if (cache$opt_lambda) {
        psi_fixed <- old_psi
        obj <- function(x) {
          ll <- cpp_tirt_item_loglik_weighted(
            cache$qt_target, cache$qt_other,
            lambda_other, psi_other, gamma_pair,
            cache$is_first, cache$Y, cache$weight, x[1L], psi_fixed)
          -ll
        }
        opt <- tryCatch(
          stats::optim(par = old_lambda, fn = obj, method = "L-BFGS-B",
                       lower = state$bounds$lambda.lower,
                       upper = state$bounds$lambda.upper,
                       control = list(maxit = state$optim.maxit)),
          error = function(e) NULL
        )
        if (!is.null(opt) && is.finite(opt$value)) {
          state$lambda[i] <- opt$par[1L]
        }

      } else if (cache$opt_psi) {
        lambda_fixed <- old_lambda
        obj <- function(x) {
          ll <- cpp_tirt_item_loglik_weighted(
            cache$qt_target, cache$qt_other,
            lambda_other, psi_other, gamma_pair,
            cache$is_first, cache$Y, cache$weight, lambda_fixed, x[1L])
          lp <- if (isTRUE(state$prior$use.prior)) {
            stats::dnorm(x[1L], mean = state$prior$psi.mu,
                         sd = state$prior$psi.sigma, log = TRUE)
          } else 0
          if (!is.finite(lp)) return(Inf)
          -ll - lp
        }
        opt <- tryCatch(
          stats::optim(par = old_psi, fn = obj, method = "L-BFGS-B",
                       lower = state$bounds$psi.lower,
                       upper = state$bounds$psi.upper,
                       control = list(maxit = state$optim.maxit)),
          error = function(e) NULL
        )
        if (!is.null(opt) && is.finite(opt$value)) {
          state$psi[i] <- opt$par[1L]
        }
      }

      max_change <- max(max_change,
                        abs(state$lambda[i] - old_lambda),
                        abs(state$psi[i]    - old_psi))
    }

    # -- Pass 2: update all gamma --
    for (cache in pair_cache) {
      if (is.null(cache)) next
      ip <- cache$ip

      # Recompute mu_no_gamma and denom (depend on lambda, psi which may have changed)
      mu_no_gamma <- state$lambda[cache$ii] * cache$qt_i -
                     state$lambda[cache$ik] * cache$qt_k
      denom <- sqrt(state$psi[cache$ii]^2 + state$psi[cache$ik]^2)

      old_gamma <- state$gamma[ip]

      obj <- function(x) {
        ll <- cpp_tirt_pair_loglik_weighted(
          mu_no_gamma, denom, cache$Y, cache$weight, x[1L])
        lp <- if (isTRUE(state$prior$use.prior)) {
          stats::dnorm(x[1L], mean = state$prior$gamma.mu,
                       sd = state$prior$gamma.sigma, log = TRUE)
        } else 0
        if (!is.finite(lp)) return(Inf)
        -ll - lp
      }

      opt <- tryCatch(
        stats::optim(par = old_gamma, fn = obj, method = "L-BFGS-B",
                     lower = state$bounds$gamma.lower,
                     upper = state$bounds$gamma.upper,
                     control = list(maxit = state$optim.maxit)),
        error = function(e) NULL
      )
      if (!is.null(opt) && is.finite(opt$value)) {
        state$gamma[ip] <- opt$par[1L]
      }

      max_change <- max(max_change, abs(state$gamma[ip] - old_gamma))
    }

    if (max_change < 1e-6) break
  }

  state
}

tirt_istem_param_vec <- function(lambda, psi, gamma, Idx.lambda.est,
                                 Idx.psi.est, Idx.gamma.est, Corr,
                                 include_corr) {
  c(lambda[Idx.lambda.est],
    psi[Idx.psi.est],
    gamma[Idx.gamma.est],
    istem_param_corr_vec(Corr, include_corr))
}

tirt_istem_param_unpack <- function(pv, pv_se, pv_rhat = NULL, state, include_corr) {
  idx <- 1L
  has_rhat <- !is.null(pv_rhat)
  lambda <- state$lambda
  lambda.se <- rep(NA_real_, length(lambda))
  lambda.rhat <- rep(NA_real_, length(lambda))
  if (length(state$Idx.lambda.est) > 0L) {
    n <- length(state$Idx.lambda.est)
    lambda[state$Idx.lambda.est] <- pv[idx:(idx + n - 1L)]
    lambda.se[state$Idx.lambda.est] <- pv_se[idx:(idx + n - 1L)]
    if (has_rhat) lambda.rhat[state$Idx.lambda.est] <- pv_rhat[idx:(idx + n - 1L)]
    idx <- idx + n
  }
  psi <- state$psi
  psi.se <- rep(NA_real_, length(psi))
  psi.rhat <- rep(NA_real_, length(psi))
  if (length(state$Idx.psi.est) > 0L) {
    n <- length(state$Idx.psi.est)
    psi[state$Idx.psi.est] <- pv[idx:(idx + n - 1L)]
    psi.se[state$Idx.psi.est] <- pv_se[idx:(idx + n - 1L)]
    if (has_rhat) psi.rhat[state$Idx.psi.est] <- pv_rhat[idx:(idx + n - 1L)]
    idx <- idx + n
  }
  gamma <- state$gamma
  gamma.se <- rep(NA_real_, length(gamma))
  gamma.rhat <- rep(NA_real_, length(gamma))
  if (length(state$Idx.gamma.est) > 0L) {
    n <- length(state$Idx.gamma.est)
    gamma[state$Idx.gamma.est] <- pv[idx:(idx + n - 1L)]
    gamma.se[state$Idx.gamma.est] <- pv_se[idx:(idx + n - 1L)]
    if (has_rhat) gamma.rhat[state$Idx.gamma.est] <- pv_rhat[idx:(idx + n - 1L)]
    idx <- idx + n
  }
  corr <- istem_unpack_corr(pv, pv_se, idx, state$D, state$Corr, include_corr)
  corr.rhat <- istem_unpack_corr_rhat(pv_rhat, idx, state$D, include_corr)
  list(
    lambda = lambda,    lambda.se = lambda.se,    lambda.rhat = lambda.rhat,
    psi = psi,          psi.se = psi.se,          psi.rhat = psi.rhat,
    gamma = gamma,      gamma.se = gamma.se,      gamma.rhat = gamma.rhat,
    Corr = corr$Corr,   Corr.se = corr$Corr.se,   Corr.rhat = corr.rhat
  )
}

tirt_istem_gamma_matrix <- function(gamma, gamma.se, pairs.matrix, I.states) {
  gamma.matrix <- matrix(NA_real_, I.states, I.states)
  gamma.matrix.se <- matrix(NA_real_, I.states, I.states)
  gamma.matrix.Rhat <- matrix(NA_real_, I.states, I.states)
  for (r in seq_len(nrow(pairs.matrix))) {
    i <- pairs.matrix[r, 1L]
    k <- pairs.matrix[r, 2L]
    gamma.matrix[i, k] <- gamma[r]
    gamma.matrix.se[i, k] <- gamma.se[r]
    gamma.matrix.Rhat[i, k] <- NA_real_
    gamma.matrix[k, i] <- -gamma[r]
    gamma.matrix.se[k, i] <- gamma.se[r]
    gamma.matrix.Rhat[k, i] <- NA_real_
  }
  list(est = gamma.matrix, se = gamma.matrix.se, Rhat = gamma.matrix.Rhat)
}
