################################# FCGDINA EM ####################################

fit.FCGDINA.EM <- function(data, Q.matrix, model, block.items, fc.type,
                           control.model, control.method, .call) {
  call <- if (is.null(.call)) match.call() else .call
  control <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  method <- fcgdina_em_method_control(control.method)
  fcgdina_em_check_method_control(method)

  prep <- fcgdina_prepare(data, Q.matrix, model, block.items, fc.type)

  bounds <- list(
    delta.lower = get_ctrl("delta.lower", -4, control.method),
    delta.upper = get_ctrl("delta.upper", 4, control.method)
  )
  if (!is.finite(bounds$delta.lower) || !is.finite(bounds$delta.upper) ||
      bounds$delta.lower >= bounds$delta.upper) {
    stop("'delta.lower' must be smaller than 'delta.upper'.", call. = FALSE)
  }

  prior.init <- rep(1.0 / prep$C, prep$C)
  set.seed(common.method$seed)

  fit <- fcgdina_em_fit_core(
    prep = prep,
    delta.init = fcgdina_initial_delta(prep),
    prior.init = prior.init,
    method = method,
    bounds = bounds,
    vis = common.method$vis
  )

  pi.est <- fit$prior

  delta.vec.est <- unlist(fit$delta, use.names = FALSE)
  delta.vec.se <- rep(NA_real_, length(delta.vec.est))
  delta.store <- NULL
  bootstrap <- NULL

  if (method$estimate.se) {
    bootstrap <- fcgdina_em_bootstrap_se(
      prep = prep,
      delta.est = fit$delta,
      prior.init = pi.est,
      method = method,
      bounds = bounds,
      vis = common.method$vis
    )
    delta.vec.se <- bootstrap$se
    if (!is.null(bootstrap$delta.mat)) {
      ok <- colSums(is.finite(bootstrap$delta.mat)) == nrow(bootstrap$delta.mat)
      if (any(ok)) {
        delta.store <- fcgdina_delta_store(bootstrap$delta.mat[, ok, drop = FALSE], prep)
      }
    }
  }

  delta.list.est <- fit$delta
  delta.list.se <- fcgdina_unpack_delta(delta.vec.se, prep)
  delta.list.Rhat <- fcgdina_unpack_delta(rep(NA_real_, length(delta.vec.est)), prep)

  prob.class <- fcgdina_prob_class(prep, delta.list.est)
  class.post.group <- fcgdina_class_posterior(prob.class, prep, pi.est,
                                              grouped = TRUE)
  class.post <- fcgdina_expand_group_matrix(class.post.group, prep)
  alpha.prob <- class.post %*% prep$alpha.patterns
  colnames(alpha.prob) <- paste0("alpha", seq_len(prep$D))
  rownames(alpha.prob) <- rownames(prep$response)
  alpha.est <- (alpha.prob > 0.5) * 1L
  storage.mode(alpha.est) <- "integer"
  colnames(class.post) <- pattern_key(prep$alpha.patterns)
  rownames(class.post) <- rownames(prep$response)

  results <- list(
    npar = prep$total_delta_len + prep$C - 1L,
    method = "EM",
    alpha = list(est = alpha.est, se = matrix(NA_real_, prep$N, prep$D),
                 Rhat = matrix(NA_real_, prep$N, prep$D),
                 prob = alpha.prob),
    class.post = class.post,
    delta = list(est = delta.list.est, se = delta.list.se,
                 Rhat = delta.list.Rhat),
    pi = pi.est,
    delta.store = delta.store,
    Q.matrix = prep$Q.matrix,
    model = prep$model,
    block.items = prep$block.items,
    patterns = prep$patterns,
    patterns.total = prep$patterns.total,
    design.matrix.list = prep$design.matrix.list,
    alpha.patterns = prep$alpha.patterns,
    fc.type = prep$fc.type,
    response = prep$response,
    stan.obj = NULL,
    MCMC.obj = NULL,
    iStEM = NULL,
    EM = list(
      iter = fit$iter,
      maxitr = method$maxitr,
      converged = fit$converged,
      tol = method$tol,
      par.tol = method$par.tol,
      logLik.trace = fit$trace,
      final.logLik = fit$logLik,
      final.ll.change = fit$ll.change,
      final.max.delta.change = fit$max.delta.change,
      bootstrap = bootstrap,
      control = method
    ),
    call = call,
    arguments = list(
      data = data,
      Q.matrix = prep$Q.matrix,
      block.items = prep$block.items,
      model = prep$model,
      fc.type = prep$fc.type,
      method = "EM",
      cores = common.method$cores,
      vis = common.method$vis,
      seed = common.method$seed,
      control.model = control,
      control.method = fit_effective_method_control(
        control.method, method = utils::modifyList(method, bounds),
        common = common.method)
    )
  )
  class(results) <- "FCGDINA"
  results$logLik <- logLik.FCGDINA(results)
  results
}

fcgdina_em_method_control <- function(control.method) {
  control.method <- fc_as_control_list(control.method, "control.method")
  estimate.se <- isTRUE(get_ctrl("estimate.se", FALSE, control.method))
  bootstrap <- as.integer(
    get_ctrl("bootstrap", if (estimate.se) 100L else 0L, control.method)
  )
  list(
    maxitr = as.integer(get_ctrl("maxitr", 500L, control.method)),
    minitr = as.integer(get_ctrl("minitr", 2L, control.method)),
    tol = get_ctrl("tol", 1e-6, control.method),
    par.tol = get_ctrl("par.tol", 1e-4, control.method),
    optim.maxit = as.integer(get_ctrl("optim.maxit", 200L, control.method)),
    estimate.se = estimate.se,
    bootstrap = bootstrap
  )
}

fcgdina_em_check_method_control <- function(method) {
  if (length(method$maxitr) != 1L || is.na(method$maxitr) ||
      method$maxitr < 1L) {
    stop("'maxitr' in 'control.method' must be a positive integer.",
         call. = FALSE)
  }
  if (length(method$minitr) != 1L || is.na(method$minitr) ||
      method$minitr < 1L || method$minitr > method$maxitr) {
    stop("'minitr' in 'control.method' must be between 1 and 'maxitr'.",
         call. = FALSE)
  }
  if (!is.numeric(method$tol) || length(method$tol) != 1L ||
      !is.finite(method$tol) || method$tol <= 0) {
    stop("'tol' in 'control.method' must be positive.", call. = FALSE)
  }
  if (!is.numeric(method$par.tol) || length(method$par.tol) != 1L ||
      !is.finite(method$par.tol) || method$par.tol <= 0) {
    stop("'par.tol' in 'control.method' must be positive.", call. = FALSE)
  }
  if (length(method$optim.maxit) != 1L || is.na(method$optim.maxit) ||
      method$optim.maxit < 1L) {
    stop("'optim.maxit' in 'control.method' must be a positive integer.",
         call. = FALSE)
  }
  if (length(method$bootstrap) != 1L || is.na(method$bootstrap) ||
      method$bootstrap < 0L) {
    stop("'bootstrap' in 'control.method' must be a non-negative integer.",
         call. = FALSE)
  }
  if (method$estimate.se && method$bootstrap < 2L) {
    stop("'bootstrap' in 'control.method' must be at least 2 when ",
         "'estimate.se = TRUE'.", call. = FALSE)
  }
  invisible(NULL)
}

fcgdina_em_fit_core <- function(prep, delta.init, prior.init, method, bounds,
                                vis = FALSE) {
  delta.current <- delta.init
  prior.current <- prior.init
  par.current <- unlist(delta.current, use.names = FALSE)
  ll.prev <- -Inf
  trace <- data.frame(
    iter = integer(0L),
    logLik = numeric(0L),
    ll.change = numeric(0L),
    max.delta.change = numeric(0L)
  )
  converged <- FALSE
  iter.done <- 0L
  ll.change <- Inf
  max.delta.change <- Inf

  if (vis) cat("\nFCGDINA EM iterations:\n")
  for (iter in seq_len(method$maxitr)) {
    # E-step: posterior class probabilities given current (delta, pi)
    prob.class <- fcgdina_prob_class(prep, delta.current)
    class.post <- fcgdina_class_posterior(prob.class, prep, prior.current,
                                          grouped = TRUE)
    class.weight <- class.post * prep$response.count

    # M-step for delta (ECM: block-wise conditional maximization)
    delta.next <- fcgdina_fit_delta(
      prep = prep,
      delta.current = delta.current,
      class.weight = class.weight,
      method = method,
      bounds = bounds
    )

    # M-step for structural parameters pi
    prior.next <- colSums(class.weight) / prep$N

    # Convergence assessment with updated parameters
    par.next <- unlist(delta.next, use.names = FALSE)
    prob.next <- fcgdina_prob_class(prep, delta.next)
    ll <- fcgdina_marginal_loglik(prob.next, prep, prior.next, grouped = TRUE)

    ll.change <- if (is.finite(ll.prev)) abs(ll - ll.prev) else Inf
    max.delta.change <- max(abs(par.next - par.current))
    trace[iter, ] <- list(iter, ll, ll.change, max.delta.change)

    if (vis) {
      istem_progress(
        "  iter %04d/%04d | logLik %12.4f | dLL %10.6f | max dpar %10.6f",
        iter, method$maxitr, ll, ll.change, max.delta.change)
    }

    delta.current <- delta.next
    prior.current <- prior.next
    par.current <- par.next
    ll.prev <- ll
    iter.done <- iter

    if (iter >= method$minitr &&
        (ll.change < method$tol || max.delta.change < method$par.tol)) {
      converged <- TRUE
      break
    }
  }
  if (vis) {
    istem_progress(
      "  iter %04d/%04d | logLik %12.4f | dLL %10.6f | max dpar %10.6f",
      iter.done, method$maxitr, ll.prev, ll.change, max.delta.change,
      newline = TRUE)
  }

  list(
    delta = delta.current,
    prior = prior.current,
    iter = iter.done,
    converged = converged,
    trace = trace,
    logLik = ll.prev,
    ll.change = ll.change,
    max.delta.change = max.delta.change
  )
}

fcgdina_em_bootstrap_se <- function(prep, delta.est, prior.init, method, bounds,
                                    vis = FALSE) {
  R <- method$bootstrap
  delta.mat <- matrix(NA_real_, nrow = prep$total_delta_len, ncol = R)
  error <- rep(NA_character_, R)

  if (vis) cat("\nFCGDINA EM bootstrap SE:\n")
  for (r in seq_len(R)) {
    idx <- sample.int(prep$N, prep$N, replace = TRUE)
    response.r <- prep$response[idx, , drop = FALSE]
    prep.r <- fcgdina_prepare(
      data = response.r,
      Q.matrix = prep$Q.matrix,
      model = prep$model,
      block.items = prep$block.items,
      fc.type = prep$fc.type
    )
    fit.r <- tryCatch(
      fcgdina_em_fit_core(
        prep = prep.r,
        delta.init = delta.est,
        prior.init = prior.init,
        method = method,
        bounds = bounds,
        vis = FALSE
      ),
      error = function(e) e
    )
    if (inherits(fit.r, "error")) {
      error[r] <- conditionMessage(fit.r)
    } else {
      delta.mat[, r] <- unlist(fit.r$delta, use.names = FALSE)
    }

    if (vis) {
      ok <- sum(colSums(is.finite(delta.mat)) == nrow(delta.mat))
      istem_progress("  bootstrap %04d/%04d | successful %04d",
                     r, R, ok)
    }
  }
  if (vis) {
    ok <- sum(colSums(is.finite(delta.mat)) == nrow(delta.mat))
    istem_progress("  bootstrap %04d/%04d | successful %04d",
                   R, R, ok, newline = TRUE)
  }

  ok <- colSums(is.finite(delta.mat)) == nrow(delta.mat)
  se <- if (sum(ok) >= 2L) {
    apply(delta.mat[, ok, drop = FALSE], 1L, stats::sd)
  } else {
    rep(NA_real_, prep$total_delta_len)
  }

  list(
    R = R,
    successful = sum(ok),
    failed = R - sum(ok),
    se = se,
    delta.mat = delta.mat,
    error = error[!is.na(error)]
  )
}
