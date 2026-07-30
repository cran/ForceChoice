#' @describeIn logLik-class FCDCM model: higher-order trait quadrature with exact
#'   marginalization over discrete attribute profiles.
#' @export
logLik.FCDCM <- function(object, theta.low = NULL, theta.up = NULL, L = NULL,
                         ...) {
  response <- object$response
  N <- nrow(response)
  B <- ncol(response)
  df <- object$npar

  grid.control <- istem_latent_grid_control(
    object$arguments$control.model, 1L, L = L,
    theta.low = theta.low, theta.up = theta.up
  )
  L <- grid.control$L
  theta.low <- grid.control$theta.lower
  theta.up <- grid.control$theta.upper

  theta.norm <- matrix(seq(theta.low, theta.up, length.out = L), ncol = 1L)
  colnames(theta.norm) <- "theta"
  pi <- stats::dnorm(theta.norm[, 1L])
  pi <- pi / sum(pi)

  delta <- object$delta$est
  support <- fcdcm_latent_support(
    theta = theta.norm,
    pi = pi,
    delta1 = delta[, "delta1"],
    delta0 = delta[, "delta0"],
    par = object$par$est,
    alpha.patterns = object$alpha.patterns,
    zeta.patterns = object$zeta.patterns,
    patterns = object$patterns
  )

  response.group <- istem_response_groups(response)
  lik <- cpp_loglik_binary(
    prob = support$prob,
    response = matrix(as.integer(response.group$response),
                      nrow = response.group$G, ncol = B),
    pi = support$pi
  )
  lik <- loglik_apply_group_count(lik, response.group$count)
  results <- c(`log Likelihood` = lik$logLik)
  results <- set_loglik_attributes(results, N, df, lik, theta.norm,
                                    support$prob, matrix(support$pi, ncol = 1L),
                                    L, theta.low, theta.up,
                                    pi.theta = matrix(pi, ncol = 1L),
                                    support.theta = support$theta,
                                    support.alpha = support$alpha,
                                    class.prob = support$class.prob,
                                    response.unique = response.group$response,
                                    response.count = response.group$count,
                                    response.group = response.group$group)
  results
}
