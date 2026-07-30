################################# FCGDINA iStEM #################################

fit.FCGDINA.iStEM <- function(data, Q.matrix, model, block.items, fc.type,
                              control.model, control.method, .call) {
  call <- if (is.null(.call)) match.call() else .call
  control <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  method <- istem_method_control(control.method)
  method <- istem_apply_grid_control(method, control, 1L)
  istem_check_method_control(method)
  control <- istem_effective_grid_model_control(control, method)

  prep <- fcgdina_prepare(data, Q.matrix, model, block.items, fc.type)

  bounds <- list(
    delta.lower = get_ctrl("delta.lower", -4, control.method),
    delta.upper = get_ctrl("delta.upper", 4, control.method)
  )
  if (!is.finite(bounds$delta.lower) || !is.finite(bounds$delta.upper) ||
      bounds$delta.lower >= bounds$delta.upper) {
    stop("'delta.lower' must be smaller than 'delta.upper'.", call. = FALSE)
  }

  delta.current <- fcgdina_initial_delta(prep)
  prior <- rep(1.0 / prep$C, prep$C)

  set.seed(common.method$seed)

  fcgdina_batch <- function(delta.current) {
    pmatrix <- matrix(NA_real_, prep$total_delta_len, method$B)
    logLik.trace <- rep(NA_real_, method$B)
    class.post <- NULL
    for (s in seq_len(method$B)) {
      prob.class <- fcgdina_prob_class(prep, delta.current)
      class.post <- fcgdina_class_posterior(prob.class, prep, prior,
                                            grouped = TRUE)
      prior <<- fcgdina_weighted_colmean(class.post, prep$response.count)
      class.weight <- fcgdina_sample_class_weight(
        class.post, count = prep$response.count)
      delta.current <- fcgdina_fit_delta(
        prep = prep, delta.current = delta.current,
        class.weight = class.weight, method = method, bounds = bounds)
      pmatrix[, s] <- unlist(delta.current, use.names = FALSE)
      prob.class <- fcgdina_prob_class(prep, delta.current)
      logLik.trace[s] <- fcgdina_marginal_loglik(
        prob.class, prep, prior, grouped = TRUE
      )
    }
    list(delta = delta.current, pmatrix = pmatrix,
         logLik.trace = logLik.trace, class.post = class.post)
  }

  run <- istem_run_batches(
    state = list(delta = delta.current, class.post = NULL),
    method = method,
    N = prep$N,
    vis = common.method$vis,
    label = "FCGDINA",
    run_batch = function(state) {
      x <- fcgdina_batch(state$delta)
      list(
        state = list(delta = x$delta, class.post = x$class.post),
        pmatrix = x$pmatrix,
        logLik.trace = x$logLik.trace
      )
    }
  )

  chain.mat <- Reduce(cbind, run$stable.plist)
  chain.mean <- rowMeans(chain.mat)
  chain.sd <- if (method$estimate.se && ncol(chain.mat) > 1L) {
    apply(chain.mat, 1L, stats::sd)
  } else {
    rep(NA_real_, length(chain.mean))
  }
  delta.list.est <- fcgdina_unpack_delta(chain.mean, prep)
  delta.list.se <- fcgdina_unpack_delta(chain.sd, prep)
  delta.list.Rhat <- fcgdina_unpack_delta(rep(NA_real_, length(chain.mean)), prep)
  delta.store <- fcgdina_delta_store(chain.mat, prep)

  prob.class <- fcgdina_prob_class(prep, delta.list.est)
  # Final E-step with converged delta and prior
  class.post.group <- fcgdina_class_posterior(prob.class, prep, prior,
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
    method = "iStEM",
    alpha = list(est = alpha.est, se = matrix(NA_real_, prep$N, prep$D),
                 Rhat = matrix(NA_real_, prep$N, prep$D),
                 prob = alpha.prob),
    class.post = class.post,
    delta = list(est = delta.list.est, se = delta.list.se,
                 Rhat = delta.list.Rhat),
    pi = prior,
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
    EM = NULL,
    iStEM = utils::modifyList(
      run$iStEM,
      list(
        accept.rate = 1,
        chain.mat = chain.mat,
        control = method
      )
    ),
    call = call,
    arguments = list(
      data = data,
      Q.matrix = prep$Q.matrix,
      block.items = prep$block.items,
      model = prep$model,
      fc.type = prep$fc.type,
      method = "iStEM",
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


fcgdina_initial_delta <- function(prep) {
  delta.current <- vector("list", prep$I.states)
  for (i in seq_len(prep$I.states)) {
    nc <- prep$delta.len[i]
    if (nc == 1L) {
      delta.current[[i]] <- 0.5
    } else if (nc == 2L) {
      delta.current[[i]] <- c(0.2, 0.6)
    } else {
      delta.current[[i]] <- c(0.2, rep(0.6 / (nc - 1L), nc - 1L))
    }
  }
  delta.current
}

fcgdina_unpack_delta <- function(delta.vec, prep) {
  delta.list <- vector("list", prep$I.states)
  for (i in seq_len(prep$I.states)) {
    idx <- prep$delta_offset[i] + seq_len(prep$delta.len[i])
    delta.list[[i]] <- as.numeric(delta.vec[idx])
  }
  delta.list
}

fcgdina_delta_store <- function(chain.mat, prep) {
  out <- vector("list", prep$I.states)
  for (i in seq_len(prep$I.states)) {
    idx <- prep$delta_offset[i] + seq_len(prep$delta.len[i])
    x <- t(chain.mat[idx, , drop = FALSE])
    colnames(x) <- paste0("delta", seq_len(ncol(x)))
    out[[i]] <- x
  }
  out
}

fcgdina_prob_class <- function(prep, delta.list) {
  delta.vec <- unlist(delta.list, use.names = FALSE)
  cpp_model_FCGDINA(
    as.matrix(prep$alpha.patterns), prep$flat_design,
    prep$design.rows, prep$design.cols,
    prep$design_offset, prep$delta_offset,
    delta.vec, prep$class.map.list,
    prep$patterns.total, prep$patterns,
    prep$C, prep$I.states)
}

fcgdina_delta_link <- function(object) {
  link <- object$delta.link
  if (is.null(link)) link <- "identity"
  match.arg(link, c("identity", "logit"))
}

fcgdina_prob_class_link <- function(prep, delta.list,
                                    link = c("identity", "logit")) {
  link <- match.arg(link)
  if (link == "identity") {
    return(fcgdina_prob_class(prep, delta.list))
  }

  prob.item <- compute_fcgdina_item_prob(
    alpha = prep$alpha.patterns,
    delta.list = delta.list,
    design.matrix.list = prep$design.matrix.list,
    Q.matrix = prep$Q.matrix,
    link = "logit"
  )
  forced_choice_from_agree(prob.item, prep$patterns.total, prep$patterns)
}

fcgdina_class_posterior <- function(prob.class, prep, prior,
                                    grouped = FALSE) {
  response <- if (grouped) prep$response.unique else prep$response
  N <- if (grouped) prep$G else prep$N
  cpp_fcgdina_class_posterior(
    prob.class, response, prep$patterns, prior, N, prep$B, prep$C)
}

fcgdina_expand_group_matrix <- function(x, prep) {
  x[prep$response.group, , drop = FALSE]
}

fcgdina_weighted_colmean <- function(x, count) {
  colSums(x * count) / sum(count)
}

fcgdina_sample_class_weight <- function(class.post, count = NULL) {
  N <- nrow(class.post)
  C <- ncol(class.post)
  cum.post <- class.post
  if (C > 1L) {
    for (j in 2L:C) {
      cum.post[, j] <- cum.post[, j] + cum.post[, j - 1L]
    }
  }
  cum.post[, C] <- 1
  cls <- max.col(cum.post >= stats::runif(N), ties.method = "first")
  weight <- matrix(0, N, C)
  weight[cbind(seq_len(N), cls)] <- if (is.null(count)) 1 else count
  weight
}

fcgdina_marginal_loglik <- function(prob.class, prep, prior,
                                    grouped = TRUE) {
  response <- if (grouped) prep$response.unique else prep$response
  count <- if (grouped) prep$response.count else rep.int(1L, nrow(response))
  log.prob <- log(pmax(prob.class, 1e-16))
  log.prior <- log(pmax(prior, 1e-16))
  ll <- 0.0

  for (n in seq_len(nrow(response))) {
    lp <- log.prior
    col.idx <- 0L
    for (b in seq_len(prep$B)) {
      y <- response[n, b]
      lp <- lp + log.prob[, col.idx + y]
      col.idx <- col.idx + prep$n.obs[b]
    }
    max.lp <- max(lp)
    ll <- ll + count[n] * (max.lp + log(sum(exp(lp - max.lp))))
  }
  ll
}

fcgdina_fit_delta <- function(prep, delta.current, class.weight,
                              method, bounds) {
  for (b in seq_len(prep$B)) {
    delta.current <- fcgdina_fit_delta_block(
      prep = prep, delta.current = delta.current,
      class.weight = class.weight, method = method,
      bounds = bounds, b = b)
  }
  delta.current
}

fcgdina_fit_delta_block <- function(prep, delta.current, class.weight,
                                    method, bounds, b) {
  items <- prep$block.items[[b]]
  delta.vec <- unlist(delta.current[items], use.names = FALSE)
  opt <- stats::optim(
    par = delta.vec,
    fn = function(par) {
      delta.block <- fcgdina_unpack_delta_block(par, prep, items)
      prob.block <- fcgdina_prob_class_block(prep, delta.block, b)
      ell <- fcgdina_weighted_loglik_block(prob.block, prep, class.weight, b)
      if (!is.finite(ell)) return(1e10)
      -ell
    },
    method = "L-BFGS-B",
    lower = rep(bounds$delta.lower, length(delta.vec)),
    upper = rep(bounds$delta.upper, length(delta.vec)),
    control = list(maxit = method$optim.maxit, fnscale = 1)
  )
  delta.current[items] <- fcgdina_unpack_delta_block(opt$par, prep, items)
  delta.current
}

fcgdina_unpack_delta_block <- function(delta.vec, prep, items) {
  delta.block <- vector("list", length(items))
  idx <- 1L
  for (j in seq_along(items)) {
    nc <- prep$delta.len[items[j]]
    delta.block[[j]] <- as.numeric(delta.vec[idx:(idx + nc - 1L)])
    idx <- idx + nc
  }
  delta.block
}

fcgdina_prob_class_block <- function(prep, delta.block, b) {
  delta.vec <- unlist(delta.block, use.names = FALSE)
  cpp_model_FCGDINA(
    as.matrix(prep$alpha.patterns),
    prep$block.flat.design[[b]],
    prep$block.design.rows[[b]],
    prep$block.design.cols[[b]],
    prep$block.design.offset[[b]],
    prep$block.delta.offset[[b]],
    delta.vec,
    prep$block.class.map.list[[b]],
    list(prep$patterns.total.local[[b]]),
    list(prep$patterns.local[[b]]),
    prep$C,
    length(prep$block.items[[b]]))
}

fcgdina_weighted_loglik_block <- function(prob.block, prep, class.weight, b) {
  N <- nrow(class.weight)
  response <- if (N == prep$G) prep$response.unique[, b] else prep$response[, b]
  cpp_fcgdina_weighted_loglik_block(
    prob.block, response, class.weight, N, prep$C)
}

# ---- Shared FCGDINA preparation -----------------------------------------------

fcgdina_prepare <- function(data, Q.matrix, model, block.items, fc.type) {
  data.matrix <- as.matrix(data)
  N <- nrow(data.matrix)
  B <- ncol(data.matrix)
  is.pattern.index.data <- is.numeric(data.matrix) || is.integer(data.matrix)

  model <- toupper(model[1L])
  if (!model %in% c("DINA", "DINO", "ACDM", "GDINA")) {
    stop("'model' must be one of: DINA, DINO, ACDM, GDINA.", call. = FALSE)
  }

  fc.type <- normalize_fc_type(fc.type, B)

  if (is.null(block.items)) {
    if (is.pattern.index.data) {
      stop("'block.items' is required when 'data' is an integer matrix of pattern indices.", call. = FALSE)
    }
    if (any(fc.type != "RANK")) {
      stop("'block.items' is required for MOLE/PICK data because partial rankings do not identify all block items.", call. = FALSE)
    }
    block.items <- get.block.items.from.data(data)
  }
  if (!is.list(block.items) || length(block.items) != B) {
    stop("'block.items' must be a list of length ", B, ".", call. = FALSE)
  }
  block.items <- lapply(block.items, as.integer)

  all.items <- unlist(block.items, use.names = FALSE)
  I.states <- length(all.items)
  if (anyDuplicated(all.items)) {
    stop("'block.items' contains duplicate items.", call. = FALSE)
  }
  if (!setequal(all.items, seq_len(I.states))) {
    stop("'block.items' must contain exactly item indices 1:", I.states, ".", call. = FALSE)
  }

  Q.matrix <- as.matrix(Q.matrix)
  if (anyNA(Q.matrix) || !all(Q.matrix %in% c(0, 1))) {
    stop("'Q.matrix' must be a 0/1 matrix.", call. = FALSE)
  }
  if (nrow(Q.matrix) != I.states) {
    stop("'Q.matrix' must have ", I.states, " rows.", call. = FALSE)
  }
  if (any(rowSums(Q.matrix) < 1L)) {
    stop("Each statement must require at least one attribute.", call. = FALSE)
  }
  D <- ncol(Q.matrix)
  storage.mode(Q.matrix) <- "numeric"

  pat <- generate_fc_permutation_patterns(block.items, fc.type)
  patterns.total <- pat$patterns.total
  patterns <- pat$patterns

  if (is.pattern.index.data) {
    response <- matrix(as.integer(data.matrix), nrow = N, ncol = B,
                       dimnames = dimnames(data.matrix))
  } else {
    response <- get.response.from.data(data, block.items, fc.type)
  }
  for (b in seq_len(B)) {
    if (anyNA(response[, b]) ||
        any(response[, b] < 1L | response[, b] > nrow(patterns[[b]]))) {
      stop("'data' contains invalid pattern indices for block ", b,
           "; valid values are 1 through ", nrow(patterns[[b]]), ".", call. = FALSE)
    }
  }
  storage.mode(response) <- "integer"

  alpha.patterns <- as.matrix(do.call(expand.grid, rep(list(0:1), D)))
  storage.mode(alpha.patterns) <- "integer"
  colnames(alpha.patterns) <- paste0("A", seq_len(D))
  C <- nrow(alpha.patterns)

  design.matrix.list <- vector("list", I.states)
  delta.len <- integer(I.states)
  design.rows <- integer(I.states)
  design.cols <- integer(I.states)
  class.map.list <- vector("list", I.states)
  class.map.mat <- matrix(NA_integer_, I.states, C)

  for (i in seq_len(I.states)) {
    att <- which(Q.matrix[i, ] > 0)
    patterns.i <- if (length(att) > 0L) {
      attributepattern(length(att))
    } else {
      matrix(0L, 1L, 0L)
    }
    X.i <- get_design_matrix_cdm(patterns.i, model)

    design.matrix.list[[i]] <- X.i
    delta.len[i] <- ncol(X.i)
    design.rows[i] <- nrow(X.i)
    design.cols[i] <- ncol(X.i)

    if (length(att) > 0L) {
      reduced.key <- pattern_key(alpha.patterns[, att, drop = FALSE])
      full.map <- match(reduced.key, pattern_key(patterns.i))
    } else {
      full.map <- rep(1L, C)
    }
    class.map.list[[i]] <- as.integer(full.map)
    class.map.mat[i, ] <- as.integer(full.map)
  }

  total_design_entries <- sum(design.rows * design.cols)
  flat_design <- numeric(total_design_entries)
  design_offset <- integer(I.states)
  delta_offset <- integer(I.states)
  off_design <- 0L
  off_delta <- 0L
  for (i in seq_len(I.states)) {
    design_offset[i] <- off_design
    delta_offset[i] <- off_delta
    X.i <- design.matrix.list[[i]]
    nr <- nrow(X.i)
    nc <- ncol(X.i)
    flat_design[(off_design + 1L):(off_design + nr * nc)] <- as.numeric(t(X.i))
    off_design <- off_design + nr * nc
    off_delta <- off_delta + nc
  }
  total_delta_len <- off_delta

  fc.type.int <- match(fc.type, c("RANK", "MOLE", "PICK"))
  n.items <- vapply(block.items, length, integer(1L))
  n.total <- vapply(patterns.total, nrow, integer(1L))
  n.obs <- vapply(patterns, nrow, integer(1L))
  fc.len <- as.integer(ifelse(fc.type == "RANK", n.items,
                              ifelse(fc.type == "MOLE", 2L, 1L)))

  response.key <- apply(response, 1L, paste, collapse = "\r")
  response.levels <- unique(response.key)
  response.group <- match(response.key, response.levels)
  response.first <- match(response.levels, response.key)
  response.unique <- response[response.first, , drop = FALSE]
  response.count <- as.integer(
    tabulate(response.group, nbins = length(response.levels))
  )
  storage.mode(response.unique) <- "integer"

  block.items.flat <- integer(sum(n.items))
  patterns.total.flat <- integer(sum(n.total * n.items))
  patterns.obs.flat <- integer(sum(n.obs * fc.len))
  patterns.total.local <- vector("list", B)
  patterns.local <- vector("list", B)
  obs.pattern.start <- c(1L, cumsum(n.obs) + 1L)
  obs.full.start <- integer(sum(n.obs) + 1L)
  obs.full.index <- integer(sum(n.total))

  item.start <- c(1L, cumsum(n.items) + 1L)
  total.start <- c(1L, cumsum(n.total * n.items) + 1L)
  obs.start <- c(1L, cumsum(n.obs * fc.len) + 1L)

  obs.full.start[1L] <- 1L
  obs.link.pos <- 1L
  obs.global <- 1L
  for (b in seq_len(B)) {
    items <- block.items[[b]]

    pat.total <- patterns.total[[b]]
    pat.total.local <- matrix(match(pat.total, items), nrow = n.total[b])
    patterns.total.local[[b]] <- pat.total.local
    total.idx <- seq.int(total.start[b], total.start[b + 1L] - 1L)
    patterns.total.flat[total.idx] <- as.integer(t(pat.total.local))

    pat.obs <- patterns[[b]]
    pat.obs.local <- matrix(match(pat.obs, items), nrow = n.obs[b])
    patterns.local[[b]] <- pat.obs.local
    obs.idx <- seq.int(obs.start[b], obs.start[b + 1L] - 1L)
    patterns.obs.flat[obs.idx] <- as.integer(t(pat.obs.local))

    for (tt in seq_len(n.obs[b])) {
      full.idx <- switch(
        fc.type[b],
        RANK = tt,
        MOLE = which(pat.total.local[, 1L] == pat.obs.local[tt, 1L] &
                       pat.total.local[, n.items[b]] == pat.obs.local[tt, 2L]),
        PICK = which(pat.total.local[, 1L] == pat.obs.local[tt, 1L])
      )
      n.link <- length(full.idx)
      obs.full.index[obs.link.pos:(obs.link.pos + n.link - 1L)] <-
        as.integer(full.idx)
      obs.link.pos <- obs.link.pos + n.link
      obs.global <- obs.global + 1L
      obs.full.start[obs.global] <- obs.link.pos
    }

    item.idx <- seq.int(item.start[b], item.start[b + 1L] - 1L)
    block.items.flat[item.idx] <- as.integer(items)
  }
  obs.full.index <- obs.full.index[seq_len(obs.link.pos - 1L)]

  block.flat.design <- vector("list", B)
  block.design.rows <- vector("list", B)
  block.design.cols <- vector("list", B)
  block.design.offset <- vector("list", B)
  block.delta.offset <- vector("list", B)
  block.class.map.list <- vector("list", B)
  for (b in seq_len(B)) {
    items <- block.items[[b]]
    block.design.rows[[b]] <- design.rows[items]
    block.design.cols[[b]] <- design.cols[items]
    block.design.offset[[b]] <- integer(length(items))
    block.delta.offset[[b]] <- integer(length(items))
    block.class.map.list[[b]] <- class.map.list[items]

    off.design <- 0L
    off.delta <- 0L
    flat.design <- numeric(sum(design.rows[items] * design.cols[items]))
    for (j in seq_along(items)) {
      i <- items[j]
      X.i <- design.matrix.list[[i]]
      nr <- nrow(X.i)
      nc <- ncol(X.i)
      n.entries <- nr * nc
      block.design.offset[[b]][j] <- off.design
      block.delta.offset[[b]][j] <- off.delta
      flat.design[(off.design + 1L):(off.design + n.entries)] <- as.numeric(t(X.i))
      off.design <- off.design + n.entries
      off.delta <- off.delta + nc
    }
    block.flat.design[[b]] <- flat.design
  }

  list(
    data = data.matrix,
    response = response,
    response.unique = response.unique,
    response.count = response.count,
    response.group = as.integer(response.group),
    G = length(response.levels),
    N = N,
    B = B,
    I.states = I.states,
    D = D,
    C = C,
    model = model,
    Q.matrix = Q.matrix,
    block.items = block.items,
    fc.type = fc.type,
    patterns = patterns,
    patterns.total = patterns.total,
    patterns.local = patterns.local,
    patterns.total.local = patterns.total.local,
    alpha.patterns = alpha.patterns,
    design.matrix.list = design.matrix.list,
    class.map.list = class.map.list,
    delta.len = delta.len,
    design.rows = design.rows,
    design.cols = design.cols,
    flat_design = flat_design,
    design_offset = design_offset,
    delta_offset = delta_offset,
    block.flat.design = block.flat.design,
    block.design.rows = block.design.rows,
    block.design.cols = block.design.cols,
    block.design.offset = block.design.offset,
    block.delta.offset = block.delta.offset,
    block.class.map.list = block.class.map.list,
    total_design_entries = total_design_entries,
    total_delta_len = total_delta_len,
    class_map_mat = class.map.mat,
    fc.type.int = fc.type.int,
    n.items = n.items,
    n.total = n.total,
    n.obs = n.obs,
    fc.len = fc.len,
    obs.pattern.start = obs.pattern.start,
    obs.full.start = obs.full.start,
    obs.full.index = obs.full.index,
    block.items.flat = block.items.flat,
    patterns.total.flat = patterns.total.flat,
    patterns.obs.flat = patterns.obs.flat,
    item.start = item.start,
    total.start = total.start,
    obs.start = obs.start
  )
}
