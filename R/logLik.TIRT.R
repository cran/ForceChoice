#' @describeIn logLik-class TIRT model: pairwise binary probit probabilities
#'   with multivariate normal quadrature over traits.
#' @export
logLik.TIRT <- function(object, theta.low = NULL, theta.up = NULL, L = NULL,
                        ...) {

  response <- object$response
  block.items <- object$arguments$block.items
  Q.matrix <- object$Q.matrix
  par <- object$par$est
  gamma.matrix <- object$gamma.matrix$est
  pairs.matrix <- object$pairs.matrix
  pairs.value  <- object$pairs.value
  fc.type      <- object$fc.type

  N.block <- length(block.items)
  N.person <- nrow(response)

  df <- object$npar

  D <- ncol(par) - 1
  grid.control <- istem_latent_grid_control(
    object$arguments$control.model, D, L = L,
    theta.low = theta.low, theta.up = theta.up
  )
  L <- grid.control$L
  theta.low <- grid.control$theta.lower
  theta.up <- grid.control$theta.upper
  I.states <- nrow(par)
  I.pairs <- ncol(response)
  lambda <- unlist(lapply(1:I.states, function(x) {
    d <- which(Q.matrix[x, ] != 0)[1]
    return(par[x, d])
  }))
  psi2 <- pmax(par[, D+1], .Machine$double.eps)

  theta.mu <- array(
    if (is.null(object$arguments$control.model$theta.mu)) rep(0, D)
    else object$arguments$control.model$theta.mu,
    dim = D
  )

  grid <- make_quadrature_grid_mvn(
    object$Corr$est, theta.mu, D, L, theta.low, theta.up
  )

  # Compute probabilities for all possible pairs
  prob <- model.TIRT(grid$theta.norm, lambda, psi2, gamma.matrix, Q.matrix, pairs.matrix)

  # Determine if we have person-specific pair data (MOLE / PICK)
  has.pairs.value <- !is.null(pairs.value) && !all(fc.type == "RANK")

  if (has.pairs.value) {
    # MOLE / PICK: use person-specific observed pairs
    I.block.vec <- sapply(block.items, length)
    block.pair.counts <- integer(N.block)
    for (b in 1:N.block) {
      k <- I.block.vec[b]
      if (fc.type[b] == "RANK") {
        block.pair.counts[b] <- choose(k, 2)
      } else if (fc.type[b] == "MOLE") {
        block.pair.counts[b] <- 2L * k - 3L
      } else if (fc.type[b] == "PICK") {
        block.pair.counts[b] <- k - 1L
      }
    }

    response.group <- istem_response_pair_groups(response, pairs.value)
    lik <- cpp_loglik_binary_person_pairs(
      prob = prob,
      response = matrix(as.integer(response.group$response),
                        nrow = response.group$G, ncol = I.pairs),
      pi = as.vector(grid$pi),
      pairs_matrix = pairs.matrix,
      pairs_value = response.group$pairs.value,
      block_pair_counts = block.pair.counts
    )
    lik <- loglik_apply_group_count(lik, response.group$count)
    results <- c(`log Likelihood` = lik$logLik)
    results <- set_loglik_attributes(results, N.person, df, lik, grid$theta.norm,
                                      prob, grid$pi, grid$L, theta.low, theta.up,
                                      response.unique = response.group$response,
                                      response.count = response.group$count,
                                      response.group = response.group$group)
    return(results)
  } else {
    response.group <- istem_response_groups(response)
    lik <- cpp_loglik_binary(
      prob = prob,
      response = matrix(as.integer(response.group$response),
                        nrow = response.group$G, ncol = I.pairs),
      pi = as.vector(grid$pi)
    )
    lik <- loglik_apply_group_count(lik, response.group$count)
    results <- c(`log Likelihood` = lik$logLik)
    results <- set_loglik_attributes(results, N.person, df, lik, grid$theta.norm,
                                      prob, grid$pi, grid$L, theta.low, theta.up,
                                      response.unique = response.group$response,
                                      response.count = response.group$count,
                                      response.group = response.group$group)
    return(results)
  }
}
