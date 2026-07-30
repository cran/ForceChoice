#' @describeIn logLik-class FCGDINA model: marginal forced-choice probabilities over
#'   all \eqn{2^D} attribute profiles.
#' @export
logLik.FCGDINA <- function(object, ...) {
  prep <- fcgdina_prepare(
    data = object$response,
    Q.matrix = object$Q.matrix,
    model = object$model,
    block.items = object$block.items,
    fc.type = object$fc.type
  )

  prob.class <- fcgdina_prob_class_link(
    prep, object$delta$est, link = fcgdina_delta_link(object))
  pi <- object$pi
  lik <- fcgdina_marginal_loglik(prob.class, prep, pi, grouped = TRUE)

  results <- c(`log Likelihood` = lik)
  attr(results, "nobs") <- prep$N
  attr(results, "df") <- object$npar
  attr(results, "pi") <- pi
  attr(results, "prob.class") <- prob.class
  class(results) <- "logLik"
  results
}
