# ===========================================================================
# rotate.MGPCM
# ===========================================================================

#' @describeIn rotate Applies a factor rotation to the posterior mean
#'   discrimination parameters of \code{MGPCM} objects. The latent trait
#'   estimates and factor correlation matrix are transformed consistently.
#'   Category intercepts and posterior uncertainty summaries are retained
#'   from the original fit without recomputation.
#' @method rotate MGPCM
#' @export
rotate.MGPCM <- function(object,
                         method,
                         vis          = TRUE,
                         target.a     = NULL,
                         target.theta = NULL,
                         pst.W        = NULL,
                         lp.p         = 1,
                         lp.gpaiter   = 5,
                         promax_m     = 4,
                         ...) {

  call <- match.call()

  D <- ncol(object$Q.matrix)
  if (is.null(D) || D < 1L) {
    stop("Cannot determine the latent dimension from object$Q.matrix.",
         call. = FALSE)
  }

  par_dimnames   <- dimnames(object$par$est)
  theta_dimnames <- dimnames(object$theta$est)
  Corr_dimnames  <- dimnames(object$Corr$est)

  mirt_object <- object
  mirt_object$par$est <- cbind(
    object$par$est[, seq_len(D), drop = FALSE],
    matrix(0, nrow(object$par$est), 3L)
  )
  mirt_object$arguments$model <- "m2pl"
  class(mirt_object) <- "MIRT"

  rotated <- rotate.MIRT(
    object       = mirt_object,
    method       = method,
    vis          = vis,
    target.a     = target.a,
    target.theta = target.theta,
    pst.W        = pst.W,
    lp.p         = lp.p,
    lp.gpaiter   = lp.gpaiter,
    promax_m     = promax_m,
    ...
  )

  out <- object
  out$par$est[, seq_len(D)] <- rotated$object$par$est[, seq_len(D), drop = FALSE]
  out$theta$est <- rotated$object$theta$est
  out$Corr$est  <- rotated$object$Corr$est

  dimnames(out$par$est)   <- par_dimnames
  dimnames(out$theta$est) <- theta_dimnames
  dimnames(out$Corr$est)  <- Corr_dimnames
  out$rotated             <- rotated$object$rotated
  out$rotation_method     <- method

  results <- list(
    object = out,
    call = call,
    arguments = list(
      object=object,
      method=method,
      vis=vis,
      target.a=target.a,
      target.theta=target.theta,
      pst.W=pst.W,
      lp.p=lp.p,
      lp.gpaiter=lp.gpaiter,
      promax_m=promax_m,
      ...
    )
  )

  class(results) <- "RotateMGPCM"
  return(results)
}
