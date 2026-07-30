################################# FCGDINA Stan #################################

#' @noRd
fit.FCGDINA.stan <- function(data, Q.matrix, model, block.items, fc.type,
                             control.model, control.method, .call) {
  call <- if (is.null(.call)) match.call() else .call
  control <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")
  common.method <- fit_common_method_control(control.method)
  stan.method <- fc_stan_method_control(control.method)

  prep <- fcgdina_prepare(data, Q.matrix, model, block.items, fc.type)

  prior <- list(
    delta.mu = istem_scalar(
      get_ctrl("delta.mu", 0.0, control),
      "delta.mu", positive = FALSE),
    delta.sigma = istem_scalar(
      get_ctrl("delta.sigma", 1.0, control),
      "delta.sigma", positive = TRUE),
    pi.alpha = istem_scalar(
      get_ctrl("pi.prior.alpha", 1.0, control),
      "pi.prior.alpha", positive = TRUE)
  )
  identification <- fcgdina_identifiable_basis(prep)

  sm <- stanmodels[["FCGDINA"]]
  if (is.null(sm)) stop("Stan model FCGDINA not found.", call. = FALSE)

  stan.data <- list(
    N = prep$N,
    G = prep$G,
    B = prep$B,
    I = prep$I.states,
    D = prep$D,
    C = prep$C,
    y_unique = prep$response.unique,
    y_count = prep$response.count,
    alpha_patterns = prep$alpha.patterns,
    total_design_entries = prep$total_design_entries,
    flat_design = prep$flat_design,
    design_rows = prep$design.rows,
    design_cols = prep$design.cols,
    design_offset = prep$design_offset,
    total_delta_len = prep$total_delta_len,
    total_delta_free = identification$rank,
    delta_basis = identification$basis,
    delta_offset = prep$delta_offset,
    class_map = prep$class_map_mat,
    n_items = prep$n.items,
    total_block_items = length(prep$block.items.flat),
    block_items = prep$block.items.flat,
    item_start = prep$item.start,
    n_total = prep$n.total,
    total_cells_total = length(prep$patterns.total.flat),
    patterns_total = prep$patterns.total.flat,
    total_start = prep$total.start,
    n_obs = prep$n.obs,
    total_obs_patterns = sum(prep$n.obs),
    obs_pattern_start = prep$obs.pattern.start,
    total_obs_full_links = length(prep$obs.full.index),
    obs_full_start = prep$obs.full.start,
    obs_full_index = prep$obs.full.index,
    delta_prior_mu = drop(crossprod(
      identification$basis, rep(prior$delta.mu, prep$total_delta_len))),
    delta_prior_sigma = prior$delta.sigma,
    pi_prior_alpha = prior$pi.alpha
  )

  mcmc_env <- setup_mcmc_env(common.method$cores)
  on.exit(restore_mcmc_env(mcmc_env), add = TRUE)
  stan.control <- build_stan_control(stan.method$algorithm, control.method)

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

  class.prob.group <- ensure_3d(MCMC.obj$class_prob_group)
  alignment <- fcgdina_align_draws(
    as.matrix(MCMC.obj$delta), as.matrix(MCMC.obj$pi),
    class.prob.group, prep, identification
  )
  MCMC.obj$delta <- alignment$delta
  MCMC.obj$pi <- alignment$pi
  MCMC.obj$class_prob_group <- alignment$class.prob
  class.prob.group <- alignment$class.prob

  # ---- Delta parameters ----
  # Stan samples delta on the logit scale.  For each stable-chain draw we
  # convert to the probability (identity) scale so that the default delta
  # output is comparable with EM / iStEM.  Logit-scale summaries are also
  # returned as delta.logit.
  delta.draws.logit <- as.matrix(MCMC.obj$delta)          # n_iter x total_delta_len

  # Logit-scale summaries (as sampled by Stan)
  delta.logit.est <- colMeans(delta.draws.logit)
  delta.logit.se <- apply(delta.draws.logit, 2L, stats::sd)

  # Probability-scale draws and summaries
  delta.draws.prob <- cpp_fcgdina_delta_logit_to_prob(
    delta.draws.logit,
    prep$flat_design,
    prep$design.rows,
    prep$design.cols,
    prep$design_offset,
    prep$delta_offset,
    prep$I.states
  )

  delta.est <- colMeans(delta.draws.prob)
  delta.se <- apply(delta.draws.prob, 2L, stats::sd)
  aligned.rhat <- fcgdina_aligned_rhat(stan.obj, prep, identification)
  delta.Rhat <- aligned.rhat$delta

  delta.list.est <- delta.list.se <- delta.list.Rhat <-
    delta.list.logit.est <- delta.list.logit.se <-
    vector("list", prep$I.states)
  for (i in seq_len(prep$I.states)) {
    idx <- prep$delta_offset[i] + seq_len(prep$delta.len[i])
    delta.list.est[[i]]       <- delta.est[idx]
    delta.list.se[[i]]        <- delta.se[idx]
    delta.list.Rhat[[i]]      <- delta.Rhat[idx]
    delta.list.logit.est[[i]] <- delta.logit.est[idx]
    delta.list.logit.se[[i]]  <- delta.logit.se[idx]
  }

  # ---- Structural parameters pi (jointly estimated by Stan) ----
  pi.draws <- as.matrix(MCMC.obj$pi)
  pi.est <- colMeans(pi.draws)
  pi.se <- apply(pi.draws, 2L, stats::sd)
  pi.Rhat <- aligned.rhat$pi
  names(pi.est) <- pattern_key(prep$alpha.patterns)

  # ---- Class posterior probabilities (already pi-weighted by Stan) ----
  class.post.group <- apply(class.prob.group, c(2, 3), mean)
  class.post <- t(class.post.group[, prep$response.group, drop = FALSE])
  log_lik.group <- as.matrix(MCMC.obj$log_lik_group)
  log_lik <- log_lik.group[, prep$response.group, drop = FALSE]

  # ---- Alpha estimates from pi-weighted posteriors ----
  alpha.prob <- class.post %*% prep$alpha.patterns
  colnames(alpha.prob) <- paste0("alpha", seq_len(prep$D))
  rownames(alpha.prob) <- rownames(prep$response)
  alpha.est <- (alpha.prob > 0.5) * 1L
  storage.mode(alpha.est) <- "integer"

  colnames(class.post) <- pattern_key(prep$alpha.patterns)
  rownames(class.post) <- rownames(prep$response)

  results <- list(
    npar = identification$rank + prep$C - 1L,
    method = "stan",
    alpha = list(est = alpha.est, se = matrix(NA_real_, prep$N, prep$D),
                 Rhat = matrix(NA_real_, prep$N, prep$D),
                 prob = alpha.prob),
    class.post = class.post,
    delta = list(est = delta.list.est, se = delta.list.se,
                 Rhat = delta.list.Rhat),
    delta.logit = list(est = delta.list.logit.est,
                       se  = delta.list.logit.se,
                       Rhat = delta.list.Rhat),
    delta.identification = identification[c(
      "rank", "nullity", "block.rank", "block.nullity", "orientation.rank",
      "convention"
    )],
    delta.link = "identity",
    delta.store = fcgdina_delta_store(t(delta.draws.prob), prep),
    pi = pi.est,
    pi.se = pi.se,
    pi.Rhat = pi.Rhat,
    stan.obj = stan.obj,
    MCMC.obj = MCMC.obj,
    log_lik = log_lik,
    EM = NULL,
    iStEM = NULL,
    Q.matrix = prep$Q.matrix,
    model = prep$model,
    block.items = prep$block.items,
    patterns = prep$patterns,
    patterns.total = prep$patterns.total,
    design.matrix.list = prep$design.matrix.list,
    alpha.patterns = prep$alpha.patterns,
    fc.type = prep$fc.type,
    response = prep$response,
    call = call,
    arguments = list(
      data = data,
      Q.matrix = prep$Q.matrix,
      block.items = prep$block.items,
      model = prep$model,
      fc.type = prep$fc.type,
      method = "stan",
      cores = common.method$cores,
      vis = common.method$vis,
      seed = common.method$seed,
      control.model = control,
      control.method = fit_effective_method_control(
        control.method, method = stan.method, common = common.method)
    )
  )
  class(results) <- "FCGDINA"
  results$logLik <- logLik.FCGDINA(results)
  results
}


#' Construct an identifiable FCGDINA coefficient basis
#'
#' Forced-choice probabilities depend only on within-block utility contrasts.
#' This helper removes every coefficient direction that leaves all such
#' contrasts unchanged. The resulting orthonormal row-space basis defines the
#' unique minimum-norm representative used by the Stan model.
#'
#' @param prep Prepared FCGDINA data returned by code{fcgdina_prepare()}.
#' @return A list containing the global basis and blockwise rank diagnostics.
#' @noRd
fcgdina_identifiable_basis <- function(prep) {
  block.basis <- vector("list", prep$B)
  block.index <- vector("list", prep$B)
  block.rank <- block.nullity <- integer(prep$B)

  for (b in seq_len(prep$B)) {
    items <- prep$block.items[[b]]
    item.len <- prep$delta.len[items]
    item.start <- cumsum(c(1L, item.len))
    n.par <- sum(item.len)
    contrast <- matrix(0, prep$C * (length(items) - 1L), n.par)
    design <- lapply(items, function(i) {
      prep$design.matrix.list[[i]][prep$class.map.list[[i]], , drop = FALSE]
    })

    for (k in 2:length(items)) {
      rows <- (k - 2L) * prep$C + seq_len(prep$C)
      contrast[rows, item.start[1L]:(item.start[2L] - 1L)] <- -design[[1L]]
      contrast[rows, item.start[k]:(item.start[k + 1L] - 1L)] <- design[[k]]
    }

    decomp <- svd(contrast, nu = 0L, nv = min(dim(contrast)))
    tol <- max(dim(contrast)) * max(decomp$d) * .Machine$double.eps
    rank <- sum(decomp$d > tol)
    if (rank < 1L) {
      stop("FCGDINA block ", b, " has no estimable utility contrast.",
           call. = FALSE)
    }

    block.basis[[b]] <- decomp$v[, seq_len(rank), drop = FALSE]
    block.index[[b]] <- unlist(lapply(items, function(i) {
      prep$delta_offset[i] + seq_len(prep$delta.len[i])
    }), use.names = FALSE)
    block.rank[b] <- rank
    block.nullity[b] <- n.par - rank
  }

  total.rank <- sum(block.rank)
  basis <- matrix(0, prep$total_delta_len, total.rank)
  col.start <- cumsum(c(1L, block.rank))
  for (b in seq_len(prep$B)) {
    cols <- col.start[b]:(col.start[b + 1L] - 1L)
    basis[block.index[[b]], cols] <- block.basis[[b]]
  }

  orientation <- matrix(0, prep$D, prep$total_delta_len)
  for (d in seq_len(prep$D)) {
    items <- which(prep$Q.matrix[, d] == 1)
    if (length(items) == 0L) {
      stop("FCGDINA attribute ", d, " is not measured by any statement.",
           call. = FALSE)
    }
    high <- prep$alpha.patterns[, d] == 1L
    for (i in items) {
      design <- prep$design.matrix.list[[i]][prep$class.map.list[[i]], ,
                                                    drop = FALSE]
      contrast <- colMeans(design[high, , drop = FALSE]) -
        colMeans(design[!high, , drop = FALSE])
      idx <- prep$delta_offset[i] + seq_len(prep$delta.len[i])
      orientation[d, idx] <- orientation[d, idx] + contrast / length(items)
    }
  }

  orientation.free <- orientation %*% basis
  orientation.svd <- svd(orientation.free, nu = 0L,
                          nv = ncol(orientation.free))
  orientation.tol <- max(dim(orientation.free)) * max(orientation.svd$d) *
    .Machine$double.eps
  orientation.rank <- sum(orientation.svd$d > orientation.tol)
  if (orientation.rank < prep$D) {
    stop("FCGDINA cannot orient all attributes: aggregate discrimination ",
         "constraints have rank ", orientation.rank, " but D = ", prep$D,
         ".", call. = FALSE)
  }

  list(
    basis = basis,
    orientation = orientation,
    rank = total.rank,
    nullity = prep$total_delta_len - total.rank,
    block.rank = block.rank,
    block.nullity = block.nullity,
    orientation.rank = orientation.rank,
    convention = paste(
      "orthonormal within-block contrast row space (minimum norm);",
      "GDINA/ACDM draws aligned to positive aggregate discrimination"
    )
  )
}


#' Align FCGDINA posterior draws to a common attribute orientation
#'
#' @param delta Draw-by-parameter logit coefficient matrix.
#' @param pi Draw-by-class structural probability matrix.
#' @param class.prob Optional draw-by-class-by-group posterior array.
#' @param prep,identification Prepared model and identification information.
#' @return Aligned draws and the applied class permutations.
#' @noRd
fcgdina_align_draws <- function(delta, pi, class.prob = NULL,
                                prep, identification) {
  delta <- as.matrix(delta)
  pi <- as.matrix(pi)
  permutations <- matrix(rep(seq_len(prep$C), each = nrow(delta)),
                         nrow(delta), prep$C)
  if (!prep$model %in% c("GDINA", "ACDM")) {
    return(list(delta = delta, pi = pi, class.prob = class.prob,
                permutations = permutations))
  }

  alpha.key <- pattern_key(prep$alpha.patterns)
  for (s in seq_len(nrow(delta))) {
    flip <- drop(identification$orientation %*% delta[s, ]) < 0
    if (!any(flip)) next

    alpha.flipped <- prep$alpha.patterns
    for (d in which(flip)) {
      alpha.flipped[, d] <- 1L - alpha.flipped[, d]
    }
    permutation <- match(pattern_key(alpha.flipped), alpha.key)
    delta.new <- delta[s, ]

    for (i in seq_len(prep$I.states)) {
      idx <- prep$delta_offset[i] + seq_len(prep$delta.len[i])
      design <- prep$design.matrix.list[[i]]
      map <- prep$class.map.list[[i]]
      eta <- drop(design[map, , drop = FALSE] %*% delta[s, idx])
      eta <- eta[permutation]
      eta.local <- vapply(seq_len(nrow(design)), function(r) {
        eta[which(map == r)[1L]]
      }, numeric(1L))
      delta.new[idx] <- qr.solve(design, eta.local)
    }

    delta[s, ] <- drop(identification$basis %*%
                          crossprod(identification$basis, delta.new))
    pi[s, ] <- pi[s, permutation]
    if (!is.null(class.prob)) {
      class.prob[s, , ] <- class.prob[s, permutation, ]
    }
    permutations[s, ] <- permutation
  }

  list(delta = delta, pi = pi, class.prob = class.prob,
       permutations = permutations)
}


#' Compute R-hat after FCGDINA attribute-orientation alignment
#'
#' @param stan.obj Fitted Stan model.
#' @param prep,identification Prepared model and identification information.
#' @return Lists of aligned delta and pi R-hat values.
#' @noRd
fcgdina_aligned_rhat <- function(stan.obj, prep, identification) {
  raw <- rstan::extract(
    stan.obj, pars = c("delta", "pi"), permuted = FALSE,
    inc_warmup = FALSE
  )
  n.iter <- dim(raw)[1L]
  n.chain <- dim(raw)[2L]
  parameter.names <- dimnames(raw)[[3L]]
  delta.idx <- grep("^delta\\[", parameter.names)
  pi.idx <- grep("^pi\\[", parameter.names)
  delta.chain <- array(NA_real_, c(n.iter, n.chain, prep$total_delta_len))
  pi.chain <- array(NA_real_, c(n.iter, n.chain, prep$C))

  for (chain in seq_len(n.chain)) {
    aligned <- fcgdina_align_draws(
      raw[, chain, delta.idx, drop = FALSE][, 1L, ],
      raw[, chain, pi.idx, drop = FALSE][, 1L, ],
      prep = prep, identification = identification
    )
    delta.chain[, chain, ] <- aligned$delta
    pi.chain[, chain, ] <- aligned$pi
  }

  if (n.chain < 2L) {
    return(list(delta = rep(NA_real_, prep$total_delta_len),
                pi = rep(NA_real_, prep$C)))
  }
  list(
    delta = rstan::monitor(delta.chain, print = FALSE)[, "Rhat"],
    pi = rstan::monitor(pi.chain, print = FALSE)[, "Rhat"]
  )
}
