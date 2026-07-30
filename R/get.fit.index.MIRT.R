#' @describeIn get.fit.index MIRT model: binary responses with Gauss--Hermite
#'   quadrature over \eqn{D}-dimensional latent space.
#' @export
get.fit.index.MIRT <- function(object, ...) {
  model <- object$arguments$model
  response <- object$arguments$response
  par <- object$par$est
  D <- object$arguments$D
  Q.matrix <- object$Q.matrix

  logLik <- object$logLik
  pi <- as.vector(attr(logLik, "pi"))
  theta.norm <- attr(logLik, "theta.norm")

  N <- nrow(response)
  I <- ncol(response)
  npar <- object$npar

  extract.free.par <- function(par, ...) {
    par.vec <- numeric(0)
    if (model != "m1pl") {
      for (i in 1:I) {
        for (j in 1:D) {
          if (Q.matrix[i, j] == 1) par.vec <- c(par.vec, par[i, j])
        }
      }
    }
    par.vec <- c(par.vec, par[, D + 1L])
    if (model %in% c("m3pl", "m4pl")) par.vec <- c(par.vec, par[, D + 2L])
    if (model == "m4pl") par.vec <- c(par.vec, par[, D + 3L])
    return(par.vec)
  }

  expand.free.par <- function(par.vec, ...) {
    par <- matrix(0, nrow = I, ncol = D + 3L)
    idx <- 1L
    if (model == "m1pl") {
      par[, 1:D] <- 1
    } else {
      for (i in 1:I) {
        for (j in 1:D) {
          if (Q.matrix[i, j] == 1) { par[i, j] <- par.vec[idx]; idx <- idx + 1L }
          else { par[i, j] <- 0.0 }
        }
      }
    }
    for (i in 1:I) { par[i, D + 1L] <- par.vec[idx]; idx <- idx + 1L }

    if (model %in% c("m1pl", "m2pl")) {
      par[, D + 2L] <- 0
    } else {
      for (i in 1:I) { par[i, D + 2L] <- par.vec[idx]; idx <- idx + 1L }
    }

    if (model %in% c("m1pl", "m2pl", "m3pl")) {
      par[, D + 3L] <- 1
    } else {
      for (i in 1:I) { par[i, D + 3L] <- par.vec[idx]; idx <- idx + 1L }
    }
    return(par)
  }

  par.vec <- extract.free.par(par, model, Q.matrix, I, D)

  loglik.fun <- function(par.vec, ...){
    par <- expand.free.par(par.vec, ...)
    prob_new <- model.MIRT(theta.norm, par)

    eps               <- .Machine$double.xmin
    log_prob          <- log(pmax(prob_new,        eps))
    log_one_minus_prob <- log(pmax(1 - prob_new,   eps))
    log_pi            <- log(pi)

    pi.length <- nrow(theta.norm)
    log_L_theta_Xi <- matrix(0, N, pi.length)
    for (p in 1:N) {
      response.p <- matrix(response[p, ], pi.length, I, byrow = TRUE)
      log_L_theta_Xi[p, ] <- rowSums(
        response.p * log_prob + (1 - response.p) * log_one_minus_prob
      ) + log_pi
    }

    logLik.value <- sum(apply(log_L_theta_Xi, 1, function(x) {
      m <- max(x)
      m + log(sum(exp(x - m)))
    }))
    return(logLik.value)
  }

  prob.fun <- function(par.vec, ...){
    par <- expand.free.par(par.vec, ...)
    prob_out <- model.MIRT(theta.norm, par)
    return(prob_out)
  }

  good.of.fit(
    par.vec    = par.vec,
    loglik.fun = loglik.fun,
    prob.fun   = prob.fun,
    response   = response,
    npar       = npar,
    pi         = pi,
    response.type = "binary",
    ...
  )
}
