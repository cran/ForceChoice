#' @describeIn get.fit.index MGPCM model: polytomous responses with M2*
#'   category-collapsed construction.
#' @export
get.fit.index.MGPCM <- function(object, ...) {
  response <- object$arguments$response
  par <- object$par$est
  D <- object$arguments$D
  Q.matrix <- object$Q.matrix
  length.poly <- object$length.poly

  logLik.obj <- object$logLik
  pi <- as.vector(attr(logLik.obj, "pi"))
  theta.norm <- attr(logLik.obj, "theta.norm")
  N <- nrow(response)
  I <- ncol(response)
  max_poly <- max(length.poly)
  npar <- object$npar

  extract.free.par <- function(par, ...) {
    par.vec <- numeric(0)
    for (i in 1:I) {
      for (j in 1:D) {
        if (Q.matrix[i, j] == 1) par.vec <- c(par.vec, par[i, j])
      }
    }
    for (i in 1:I) {
      Ki <- length.poly[i]
      if (Ki > 1) {
        cols <- (D + 2):(D + Ki)
        par.vec <- c(par.vec, par[i, cols])
      }
    }
    return(par.vec)
  }

  expand.free.par <- function(par.vec, ...) {
    par <- matrix(NA_real_, nrow = I, ncol = D + max_poly)
    idx <- 1L
    for (i in 1:I) {
      for (j in 1:D) {
        if (Q.matrix[i, j] == 1) { par[i, j] <- par.vec[idx]; idx <- idx + 1L }
        else { par[i, j] <- 0.0 }
      }
    }
    par[, D + 1] <- 0
    for (i in 1:I) {
      Ki <- length.poly[i]
      if (Ki > 1) {
        for (k in 1:(Ki - 1)) {
          par[i, D + 1 + k] <- par.vec[idx]
          idx <- idx + 1L
        }
      }
    }
    return(par)
  }

  par.vec <- extract.free.par(par)

  response.unfold <- matrix(0, N, sum(length.poly))
  idx <- 0
  for (i in 1:I) {
    response.temp <- matrix(0, N, length.poly[i])
    for (p in 1:N) {
      response.temp[p, response[p, i] + 1] <- 1
    }
    response.unfold[, (idx + 1):(idx + length.poly[i])] <- response.temp
    idx <- idx + length.poly[i]
  }

  loglik.fun <- function(par.vec, ...) {
    par <- expand.free.par(par.vec, ...)
    prob_new <- model.MGPCM(theta.norm, par)

    eps        <- .Machine$double.xmin
    log_prob   <- log(pmax(prob_new, eps))
    log_pi     <- log(pi)

    pi.length <- nrow(theta.norm)
    log_L_theta_Xi <- matrix(0, N, pi.length)
    for (p in 1:N) {
      response.p <- matrix(response.unfold[p, ], pi.length, ncol(prob_new), byrow = TRUE)
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
    model.MGPCM(theta.norm, par)
  }

  good.of.fit(
    par.vec    = par.vec,
    loglik.fun = loglik.fun,
    prob.fun   = prob.fun,
    response   = response,
    npar       = npar,
    pi         = pi,
    response.type = "polytomous",
    response.K = length.poly,
    ...
  )
}
