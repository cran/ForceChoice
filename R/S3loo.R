# ---- Extract pointwise log-likelihood -----------------------------------------

#' Extract pointwise log-likelihood matrix
#'
#' For Stan-fitted models, returns the posterior draws of pointwise
#' log-likelihood (an \eqn{S \times N_{\text{obs}}} matrix for most models,
#' or \eqn{S \times N} for FCDCM and FCGDINA).  For EM and iStEM fits,
#' returns \code{NULL}.
#'
#' The returned matrix can be passed directly to \code{\link[loo]{loo}} or
#' \code{\link[loo]{waic}} for model comparison.
#'
#' @param object A fitted ForceChoice model object.
#' @param ...    Additional arguments (currently ignored).
#' @return An \eqn{S \times N} matrix of pointwise log-likelihood draws
#'   (Stan), or \code{NULL} (EM / iStEM).
#'
#' @examples
#' \donttest{
#' # stan code, long time
#' sim <- sim.data.MIRT(N = 20, I = 6, D = 2, model = "2PL")
#' fit <- fit.MIRT(sim$response, model = "2PL", D = 2, method = "stan")
#' log_lik <- get.log_lik(fit)
#' loo::loo(log_lik)
#' loo::waic(log_lik)
#' }
#' @export
get.log_lik <- function(object, ...) {
  UseMethod("get.log_lik")
}

#' @describeIn get.log_lik Default method; returns \code{object$log_lik}
#'   when present and \code{NULL} otherwise.
#' @export
get.log_lik.default <- function(object, ...) {
  if (!is.null(object$log_lik)) {
    return(object$log_lik)
  }
  NULL
}
