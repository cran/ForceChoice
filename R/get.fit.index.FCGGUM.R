#' @describeIn get.fit.index FCGGUM model: forced-choice unfolding with nominal
#'   binary expansion; within-block pairs excluded.
#' @export
get.fit.index.FCGGUM <- function(object, ...) {
  response <- object$response
  par <- object$par$est
  D <- object$arguments$D
  Q.matrix <- object$Q.matrix
  block.items <- object$block.items
  patterns <- object$patterns
  patterns.total <- object$patterns.total

  logLik.obj <- object$logLik
  pi <- as.vector(attr(logLik.obj, "pi"))
  theta.norm <- attr(logLik.obj, "theta.norm")
  N <- nrow(response)
  B <- length(block.items)
  I <- length(unlist(block.items))
  length.poly <- object$length.poly
  max_poly <- max(length.poly)
  npar <- object$npar

  n_cat <- sapply(patterns, nrow)
  n_bin <- sum(n_cat - 1L)
  response.bin <- matrix(0L, N, n_bin)
  indicator.block <- integer(n_bin)
  idx.resp <- 0L
  for (b in 1:B) {
    Kb <- n_cat[b]
    if (Kb > 1L) {
      for (j in 1:(Kb - 1L)) {
        response.bin[, idx.resp + j] <- as.integer(response[, b] == j)
        indicator.block[idx.resp + j] <- b
      }
      idx.resp <- idx.resp + (Kb - 1L)
    }
  }

  extract.free.par <- function(par, ...) {
    par.vec <- numeric(0)
    for (i in 1:I) {
      for (j in 1:D) {
        if (abs(Q.matrix[i, j]) == 1) par.vec <- c(par.vec, par[i, j])
      }
    }
    for (i in 1:I) {
      for (j in 1:D) {
        if (abs(Q.matrix[i, j]) == 1) par.vec <- c(par.vec, par[i, D + j])
      }
    }
    for (i in 1:I) {
      Ki <- length.poly[i]
      if (Ki > 1) {
        cols <- (D + D + 2):(D + D + Ki)
        par.vec <- c(par.vec, par[i, cols])
      }
    }
    return(par.vec)
  }

  expand.free.par <- function(par.vec, ...) {
    par <- matrix(NA_real_, nrow = I, ncol = D + D + max_poly)
    idx <- 1L
    for (i in 1:I) {
      for (j in 1:D) {
        if (abs(Q.matrix[i, j]) == 1) { par[i, j] <- par.vec[idx]; idx <- idx + 1L }
        else { par[i, j] <- 0.0 }
      }
    }
    for (i in 1:I) {
      for (j in 1:D) {
        if (abs(Q.matrix[i, j]) == 1) {
          par[i, D + j] <- par.vec[idx]
          idx <- idx + 1L
        } else {
          par[i, D + j] <- 0.0
        }
      }
    }
    par[, D + D + 1] <- 0
    for (i in 1:I) {
      Ki <- length.poly[i]
      if (Ki > 1) {
        for (k in 1:(Ki - 1)) {
          par[i, D + D + 1 + k] <- par.vec[idx]
          idx <- idx + 1L
        }
      }
    }
    return(par)
  }

  par.vec <- extract.free.par(par)

  response.full <- matrix(0, N, sum(n_cat))
  idx.full <- 0L
  for (b in 1:B) {
    Kb <- n_cat[b]
    for (p in 1:N) {
      response.full[p, idx.full + response[p, b]] <- 1
    }
    idx.full <- idx.full + Kb
  }

  loglik.fun <- function(par.vec, ...) {
    par <- expand.free.par(par.vec, ...)
    prob_new <- model.FCGGUM(theta.norm, par, patterns.total, patterns)

    eps        <- .Machine$double.xmin
    log_prob   <- log(pmax(prob_new, eps))
    log_pi     <- log(pi)

    pi.length <- nrow(theta.norm)
    log_L_theta_Xi <- matrix(0, N, pi.length)
    for (p in 1:N) {
      response.p <- matrix(response.full[p, ], pi.length, ncol(prob_new), byrow = TRUE)
      log_L_theta_Xi[p, ] <- rowSums(response.p * log_prob) + log_pi
    }

    logLik.value <- sum(apply(log_L_theta_Xi, 1, function(x) {
      m <- max(x)
      m + log(sum(exp(x - m)))
    }))
    return(logLik.value)
  }

  prob.fun <- function(par.vec, ...) {
    par <- expand.free.par(par.vec, ...)
    prob_full <- model.FCGGUM(theta.norm, par, patterns.total, patterns)
    prob_out <- matrix(0, nrow(prob_full), n_bin)
    idx.full <- 0L
    idx.bin  <- 0L
    for (b in 1:B) {
      Kb <- n_cat[b]
      if (Kb > 1L) {
        prob_out[, (idx.bin + 1L):(idx.bin + Kb - 1L)] <-
          prob_full[, (idx.full + 1L):(idx.full + Kb - 1L)]
        idx.bin  <- idx.bin + (Kb - 1L)
      }
      idx.full <- idx.full + Kb
    }
    return(prob_out)
  }

  good.of.fit(
    par.vec          = par.vec,
    loglik.fun       = loglik.fun,
    prob.fun         = prob.fun,
    response         = response.bin,
    npar             = npar,
    pi               = pi,
    response.type    = "nominal",
    nominal.groups   = indicator.block,
    ...
  )
}
