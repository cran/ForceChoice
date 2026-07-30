#' @describeIn get.fit.index TIRT model: pairwise binary responses with
#'   within-block pairs excluded via \code{bivariate.groups}.
#' @export
get.fit.index.TIRT <- function(object, ...) {
  block.items <- object$arguments$block.items
  Q.matrix <- object$Q.matrix
  response <- object$response
  par <- object$par$est
  D <- ncol(par) - 1L

  gamma.matrix <- object$gamma.matrix$est
  pairs.matrix <- object$pairs.matrix
  pairs.value  <- object$pairs.value
  fc.type      <- object$fc.type

  logLik <- object$logLik
  pi <- as.vector(attr(logLik, "pi"))
  theta.norm <- attr(logLik, "theta.norm")
  N.person <- nrow(response)
  I.pairs <- ncol(response)
  I.states <- nrow(par)
  npar <- object$npar
  N.block <- length(block.items)

  I.block <- length(block.items[[1]])
  is_case_a <- (I.block == 2 && D > 2)
  is_case_b <- (D == 2 && I.block == 2)

  if (is_case_b) {
    idx_lambda_fix <- block.items[[1]]
  } else {
    idx_lambda_fix <- integer(0)
  }

  if (is_case_a || is_case_b) {
    idx_psi_fix <- 1:I.states
  } else {
    idx_psi_fix <- sapply(block.items, function(x) x[length(x)])
  }

  idx_psi_est <- setdiff(1:I.states, idx_psi_fix)

  pairs.char <- apply(pairs.matrix, 1, function(x){
    return(paste0(x[1], x[2]))
  })

  I.block.vec <- sapply(block.items, length)
  pairs.per.block.all <- choose(I.block.vec, 2L)
  block.pair.counts <- integer(N.block)
  for (b in seq_len(N.block)) {
    k <- I.block.vec[b]
    if (fc.type[b] == "RANK") {
      block.pair.counts[b] <- choose(k, 2L)
    } else if (fc.type[b] == "MOLE") {
      block.pair.counts[b] <- 2L * k - 3L
    } else if (fc.type[b] == "PICK") {
      block.pair.counts[b] <- k - 1L
    }
  }

  has.pairs.value <- !is.null(pairs.value) && !all(fc.type == "RANK")

  if (has.pairs.value) {
    expanded <- cpp_tirt_expand_response_full(
      response = matrix(as.integer(response), nrow = N.person, ncol = ncol(response)),
      pairs_matrix = pairs.matrix,
      pairs_value = pairs.value,
      block_pair_counts = block.pair.counts
    )
    response.full <- expanded$response.full
    response.use <- response.full
    idx_gamma_est <- expanded$idx.gamma.est
  } else {
    response.use <- response
    idx_gamma_est <- seq_len(nrow(pairs.matrix))
  }
  idx_gamma_fix <- setdiff(seq_len(nrow(pairs.matrix)), idx_gamma_est)

  if (has.pairs.value) {
    cb.resp <- rep(seq_len(N.block), times = pairs.per.block.all)
  } else {
    cb.resp <- rep(seq_len(N.block), times = block.pair.counts)
  }

  extract.free.par <- function(par, Q.matrix, idx_lambda_fix, idx_psi_est,
                               gamma.matrix, pairs.matrix, D, pairs.char) {
    par.vec <- numeric(0)
    for(i in 1:I.states){
      if(! i %in% idx_lambda_fix){
        d <- which(Q.matrix[i, ] != 0)[1]
        par.vec <- c(par.vec, par[i, d])
      }
    }
    for (i in idx_psi_est) {
      par.vec <- c(par.vec, par[i, D + 1L])
    }

    if (length(idx_gamma_est) > 0) {
      pairs.char.est <- pairs.char[idx_gamma_est]
      for (i in seq_along(pairs.char.est)) {
        pos <- match(pairs.char.est[i], pairs.char)
        par.vec <- c(par.vec, gamma.matrix[pairs.matrix[pos, 1], pairs.matrix[pos, 2]])
      }
    }

    return(par.vec)
  }

  expand.free.par <- function(par.vec, ...) {
    lambda <- rep(NA_real_, I.states)
    lambda[idx_lambda_fix] <- 0.80
    if(is_case_a || is_case_b){
      psi2 <- rep(0.50, I.states)
    }else{
      psi2 <- rep(1.00, I.states)
    }
    gamma.matrix <- matrix(NA_real_, nrow = I.states, ncol = I.states)

    idx <- 1L
    for (i in 1:I.states) {
      if (! i %in% idx_lambda_fix) {
        lambda[i] <- par.vec[idx]
        idx <- idx + 1L
      }
    }

    for (i in idx_psi_est) {
      psi2[i] <- par.vec[idx]
      idx <- idx + 1L
    }

    for (fix_idx in idx_gamma_fix) {
      r <- pairs.matrix[fix_idx, 1]
      c <- pairs.matrix[fix_idx, 2]
      gamma.matrix[r, c] <- 0
      gamma.matrix[c, r] <- 0
    }
    if (length(idx_gamma_est) > 0) {
      for (i in seq_along(idx_gamma_est)) {
        r <- pairs.matrix[idx_gamma_est[i], 1]
        c <- pairs.matrix[idx_gamma_est[i], 2]
        gamma.matrix[r, c] <- par.vec[idx + i - 1]
        gamma.matrix[c, r] <- -par.vec[idx + i - 1]
      }
    }
    psi2 <- pmax(psi2, .Machine$double.eps)
    return(list(lambda = lambda, psi2 = psi2, gamma.matrix = gamma.matrix))
  }

  par.vec <- extract.free.par(par, Q.matrix, idx_lambda_fix, idx_psi_est,
                              gamma.matrix, pairs.matrix, D, pairs.char)

  loglik.fun <- function(par.vec, ...) {
    mats <- expand.free.par(par.vec, ...)
    prob_new <- model.TIRT(theta.norm, mats$lambda, mats$psi2,
                           mats$gamma.matrix, Q.matrix, pairs.matrix)

    if (has.pairs.value) {
      lik <- cpp_loglik_binary_person_pairs(
        prob = prob_new,
        response = matrix(as.integer(response), nrow = N.person, ncol = ncol(response)),
        pi = pi,
        pairs_matrix = pairs.matrix,
        pairs_value = pairs.value,
        block_pair_counts = block.pair.counts
      )
    } else {
      lik <- cpp_loglik_binary(
        prob = prob_new,
        response = matrix(as.integer(response), nrow = N.person, ncol = I.pairs),
        pi = pi
      )
    }

    lik$logLik
  }

  prob.fun <- function(par.vec, ...) {
    mats <- expand.free.par(par.vec, ...)
    model.TIRT(theta.norm, mats$lambda, mats$psi2,
               mats$gamma.matrix, Q.matrix, pairs.matrix)
  }

  good.of.fit(
    par.vec     = par.vec,
    loglik.fun  = loglik.fun,
    prob.fun    = prob.fun,
    response    = response.use,
    npar        = npar,
    pi          = pi,
    response.type = "binary",
    bivariate.groups = cb.resp,
    ...
  )
}
