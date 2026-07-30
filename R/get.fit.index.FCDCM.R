#' @describeIn get.fit.index FCDCM model: binary block responses with
#'   quadrature over higher-order \eqn{\theta} and exact marginalization
#'   over \eqn{2^D} attribute patterns.
#' @export
get.fit.index.FCDCM <- function(object, ...) {
  response <- object$response
  D <- object$arguments$D
  B <- ncol(response)
  npar <- object$npar

  logLik.obj <- object$logLik
  pi <- as.vector(attr(logLik.obj, "pi"))
  theta.norm <- attr(logLik.obj, "theta.norm")
  pi.theta <- as.vector(attr(logLik.obj, "pi.theta"))
  if (length(pi.theta) == 0L) {
    pi.theta <- stats::dnorm(theta.norm[, 1L])
    pi.theta <- pi.theta / sum(pi.theta)
  }

  par.vec <- fcdcm_free_par_vec(
    delta1 = object$delta$est[, "delta1"],
    delta0 = object$delta$est[, "delta0"],
    par = object$par$est
  )

  loglik.fun <- function(par.vec, ...) {
    pars <- fcdcm_expand_free_par(
      par.vec, D = D, B = B,
      dim.names = rownames(object$delta$est),
      block.names = rownames(object$par$est),
      par.names = colnames(object$par$est)
    )
    support <- fcdcm_latent_support(
      theta = theta.norm,
      pi = pi.theta,
      delta1 = pars$delta1,
      delta0 = pars$delta0,
      par = pars$par,
      alpha.patterns = object$alpha.patterns,
      zeta.patterns = object$zeta.patterns,
      patterns = object$patterns
    )
    cpp_loglik_binary(
      prob = support$prob,
      response = matrix(as.integer(response), nrow = nrow(response), ncol = B),
      pi = support$pi
    )$logLik
  }

  prob.fun <- function(par.vec, ...) {
    pars <- fcdcm_expand_free_par(
      par.vec, D = D, B = B,
      dim.names = rownames(object$delta$est),
      block.names = rownames(object$par$est),
      par.names = colnames(object$par$est)
    )
    fcdcm_latent_support(
      theta = theta.norm,
      pi = pi.theta,
      delta1 = pars$delta1,
      delta0 = pars$delta0,
      par = pars$par,
      alpha.patterns = object$alpha.patterns,
      zeta.patterns = object$zeta.patterns,
      patterns = object$patterns
    )$prob
  }

  pi.fun <- function(par.vec, ...) {
    pars <- fcdcm_expand_free_par(
      par.vec, D = D, B = B,
      dim.names = rownames(object$delta$est),
      block.names = rownames(object$par$est),
      par.names = colnames(object$par$est)
    )
    fcdcm_latent_support(
      theta = theta.norm,
      pi = pi.theta,
      delta1 = pars$delta1,
      delta0 = pars$delta0,
      par = pars$par,
      alpha.patterns = object$alpha.patterns,
      zeta.patterns = object$zeta.patterns,
      patterns = object$patterns
    )$pi
  }

  good.of.fit(
    par.vec = par.vec,
    loglik.fun = loglik.fun,
    prob.fun = prob.fun,
    response = response,
    npar = npar,
    pi = pi,
    pi.fun = pi.fun,
    response.type = "binary",
    ...
  )
}
