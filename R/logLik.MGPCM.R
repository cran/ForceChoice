#' @describeIn logLik-class MGPCM model: polytomous responses integrated over the
#'   fitted multivariate normal latent distribution.
#' @export
logLik.MGPCM <- function(object, theta.low = NULL, theta.up = NULL, L = NULL,
                         ...) {

  par <- object$par$est
  response <- object$arguments$response
  length.poly <- object$length.poly

  N <- nrow(response)
  I <- ncol(response)

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

  prob <- model.MGPCM(grid$theta.norm, par)
  response.group <- istem_response_groups(response)
  lik <- cpp_loglik_indexed(
    prob = prob,
    response = matrix(as.integer(response.group$response),
                      nrow = response.group$G, ncol = I),
    pi = as.vector(grid$pi),
    block_sizes = as.integer(length.poly),
    response_base = 0L
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
