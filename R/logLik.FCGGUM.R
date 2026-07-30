#' @describeIn logLik-class FCGGUM model: forced-choice unfolding block
#'   probabilities with multivariate normal quadrature.
#' @export
logLik.FCGGUM <- function(object, theta.low = NULL, theta.up = NULL, L = NULL,
                          ...) {

  par <- object$par$est
  response <- object$response
  block.items <- object$block.items
  patterns <- object$patterns
  patterns.total <- object$patterns.total

  N <- nrow(response)
  B <- length(block.items)

  df <- object$npar
  D <- object$arguments$D
  grid.control <- istem_latent_grid_control(
    object$arguments$control.model, D, L = L,
    theta.low = theta.low, theta.up = theta.up
  )
  L <- grid.control$L
  theta.low <- grid.control$theta.lower
  theta.up <- grid.control$theta.upper

  theta.mu <- array(
    if (is.null(object$arguments$control.model$theta.mu)) rep(0, D)
    else object$arguments$control.model$theta.mu,
    dim = D
  )

  grid <- make_quadrature_grid_mvn(
    object$Corr$est, theta.mu, D, L, theta.low, theta.up
  )

  prob <- model.FCGGUM(grid$theta.norm, par, patterns.total, patterns)
  response.group <- istem_response_groups(response)
  lik <- cpp_loglik_indexed(
    prob = prob,
    response = matrix(as.integer(response.group$response),
                      nrow = response.group$G, ncol = B),
    pi = as.vector(grid$pi),
    block_sizes = as.integer(vapply(patterns, nrow, integer(1L))),
    response_base = 1L
  )
  lik <- loglik_apply_group_count(lik, response.group$count)
  results <- c(`log Likelihood` = lik$logLik)
  results <- set_loglik_attributes(results, N, df, lik, grid$theta.norm,
                                    prob, grid$pi, grid$L, theta.low, theta.up,
                                    response.unique = response.group$response,
                                    response.count = response.group$count,
                                    response.group = response.group$group)
  return(results)
}
