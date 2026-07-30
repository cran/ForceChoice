#' @describeIn get.fit.index FCMIRT model: forced-choice ranking data expanded
#'   to nominal binary indicators; within-block pairs excluded from bivariate
#'   moments.
#' @export
get.fit.index.FCMIRT <- function(object, ...) {
  model <- object$arguments$model
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
  I <- length(unlist(block.items, use.names = FALSE))
  npar <- object$npar

  n_cat <- vapply(patterns, nrow, integer(1L))
  n_bin <- sum(n_cat - 1L)
  response.bin <- matrix(0L, N, n_bin)
  indicator.block <- integer(n_bin)

  idx.resp <- 0L
  for (b in seq_len(B)) {
    Kb <- n_cat[b]
    if (Kb > 1L) {
      for (j in seq_len(Kb - 1L)) {
        response.bin[, idx.resp + j] <- as.integer(response[, b] == j)
        indicator.block[idx.resp + j] <- b
      }
      idx.resp <- idx.resp + (Kb - 1L)
    }
  }

  extract.free.par <- function(par, ...) {
    par.vec <- numeric(0)
    if (model != "m1pl") {
      for (i in seq_len(I)) {
        for (j in seq_len(D)) {
          if (Q.matrix[i, j] == 1) par.vec <- c(par.vec, par[i, j])
        }
      }
    }
    for (items in block.items) {
      par.vec <- c(par.vec, par[items[-length(items)], D + 1L])
    }
    if (model %in% c("m3pl", "m4pl")) par.vec <- c(par.vec, par[, D + 2L])
    if (model == "m4pl") par.vec <- c(par.vec, par[, D + 3L])
    return(par.vec)
  }

  expand.free.par <- function(par.vec, ...) {
    par <- matrix(0, nrow = I, ncol = D + 3L)
    idx <- 1L
    if (model == "m1pl") {
      par[, seq_len(D)] <- Q.matrix
    } else {
      for (i in seq_len(I)) {
        for (j in seq_len(D)) {
          if (Q.matrix[i, j] == 1) {
            par[i, j] <- par.vec[idx]
            idx <- idx + 1L
          } else {
            par[i, j] <- 0.0
          }
        }
      }
    }

    for (items in block.items) {
      K <- length(items)
      b_free <- par.vec[idx:(idx + K - 2L)]
      par[items[-K], D + 1L] <- b_free
      par[items[K], D + 1L] <- -sum(b_free)
      idx <- idx + K - 1L
    }

    if (model %in% c("m1pl", "m2pl")) {
      par[, D + 2L] <- 0
    } else {
      for (i in seq_len(I)) {
        par[i, D + 2L] <- par.vec[idx]
        idx <- idx + 1L
      }
    }

    if (model %in% c("m1pl", "m2pl", "m3pl")) {
      par[, D + 3L] <- 1
    } else {
      for (i in seq_len(I)) {
        par[i, D + 3L] <- par.vec[idx]
        idx <- idx + 1L
      }
    }
    return(par)
  }

  par.vec <- extract.free.par(par)

  response.full <- matrix(0, N, sum(n_cat))
  idx.full <- 0L
  for (b in seq_len(B)) {
    Kb <- n_cat[b]
    for (p in seq_len(N)) {
      response.full[p, idx.full + response[p, b]] <- 1
    }
    idx.full <- idx.full + Kb
  }

  loglik.fun <- function(par.vec, ...) {
    par <- expand.free.par(par.vec, ...)
    prob_new <- model.FCMIRT(theta.norm, par, patterns.total, patterns)

    log_prob <- log(pmax(prob_new, .Machine$double.xmin))
    log_pi <- log(pi)

    pi.length <- nrow(theta.norm)
    log_L_theta_Xi <- matrix(0, N, pi.length)
    for (p in seq_len(N)) {
      response.p <- matrix(response.full[p, ], pi.length, ncol(prob_new), byrow = TRUE)
      log_L_theta_Xi[p, ] <- rowSums(response.p * log_prob) + log_pi
    }

    sum(apply(log_L_theta_Xi, 1, function(x) {
      m <- max(x)
      m + log(sum(exp(x - m)))
    }))
  }

  prob.fun <- function(par.vec, ...) {
    par <- expand.free.par(par.vec, ...)
    prob_full <- model.FCMIRT(theta.norm, par, patterns.total, patterns)
    prob_out <- matrix(0, nrow(prob_full), n_bin)

    idx.full <- 0L
    idx.bin <- 0L
    for (b in seq_len(B)) {
      Kb <- n_cat[b]
      if (Kb > 1L) {
        prob_out[, (idx.bin + 1L):(idx.bin + Kb - 1L)] <-
          prob_full[, (idx.full + 1L):(idx.full + Kb - 1L)]
        idx.bin <- idx.bin + (Kb - 1L)
      }
      idx.full <- idx.full + Kb
    }
    return(prob_out)
  }

  good.of.fit(
    par.vec = par.vec,
    loglik.fun = loglik.fun,
    prob.fun = prob.fun,
    response = response.bin,
    npar = npar,
    pi = pi,
    response.type = "nominal",
    nominal.groups = indicator.block,
    ...
  )
}
