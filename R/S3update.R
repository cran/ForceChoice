#' @title S3 Methods: update
#'
#' @description
#' The \code{update} function provides a unified and convenient interface to
#' re-fit \pkg{ForceChoice} models with modified parameter settings while
#' preserving all other original configurations. It allows users to change
#' the estimation method, adjust control parameters, increase dimensionality,
#' or switch model specifications without re-specifying the entire call.
#'
#' @param object An object of one of the following classes:
#'   \itemize{
#'     \item \code{"MIRT"} --- Multidimensional IRT model.
#'     \item \code{"MGPCM"} --- Multidimensional Generalized Partial Credit Model.
#'     \item \code{"MGGUM"} --- Multidimensional Generalized Graded Unfolding Model.
#'     \item \code{"FCMIRT"} --- Forced-Choice Multidimensional IRT model.
#'     \item \code{"FCDCM"} --- Forced-Choice Diagnostic Classification Model.
#'     \item \code{"FCGDINA"} --- Forced-Choice GDINA model.
#'     \item \code{"FCGGUM"} --- Forced-Choice Generalized Graded Unfolding Model.
#'     \item \code{"TIRT"} --- Thurstonian IRT for Forced-Choice.
#'   }
#' @param ... Additional named arguments passed to override or extend the
#'   original call. Valid arguments depend on the class of \code{object} and
#'   correspond to the formal parameters of the matching \code{fit.*()}
#'   function. Common overrides include \code{method}, \code{control.model},
#'   \code{control.method}, and legacy \code{control}.
#'
#' @return An object of the same class as \code{object}, re-fitted using the
#'   original arguments updated with any provided in \code{...}. All
#'   unchanged parameters are preserved from the original call.
#'
#' @details
#' Internally, each method extracts the stored \code{arguments} list from
#' the input \code{object}, merges it with user-provided \code{...} using
#' \code{\link[utils]{modifyList}}, explicitly merges control lists from the
#' previous fit and the update call, filters to only the formal parameters of
#' the target \code{fit.*()} function, then re-invokes the fitting function
#' via \code{do.call}.
#'
#' This ensures that:
#' \itemize{
#'   \item Only explicitly overridden parameters are changed.
#'   \item Default values from the original call remain intact.
#'   \item Complex nested structures (e.g., \code{control.model},
#'         \code{control.method}, and \code{control} lists) can be partially
#'         updated.
#'   \item Internal arguments (e.g., \code{cores}, \code{vis}, \code{seed})
#'         are passed through to the method controls.
#' }
#'
#' @examples
#' sim <- sim.data.MIRT(N = 20, I = 6, D = 2, model = "m2pl")
#' fit <- fit.MIRT(sim$response, model = "m2pl", D = 2, method = "iStEM",
#'                 control.method = list(
#'                   vis = FALSE, seed = 123,
#'                   M = 2, B = 2, burnin.maxitr = 2,
#'                   maxitr = 3, eps1 = 10, eps2 = 10,
#'                   estimate.se = FALSE))
#'
#' \donttest{
#' # stan code, long time
#' # Switch to Stan estimation with more chains
#' fit_stan <- update(fit, method = "stan",
#'                    control.method = list(chains = 1, iter = 200))
#' }
#'
#' # Change model to 3PL
#' fit_3pl <- update(fit, model = "m3pl")
#'
#' # Adjust prior hyperparameters
#' fit_prior <- update(fit, control.model = list(a.mu = 0.5, a.sigma = 0.5))
#'
#' @name update
#' @importFrom stats update
NULL

# ---- Internal helper ----

.fc_update_args <- function(arguments, fun, ...) {
  if (is.null(arguments$data) && !is.null(arguments$response)) {
    arguments$data <- arguments$response
  }
  arguments$response <- NULL

  dots <- list(...)
  control.names <- c("control.model", "control.method", "control")
  for (nm in intersect(control.names, union(names(arguments), names(dots)))) {
    old <- fc_as_control_list(arguments[[nm]], nm)
    new <- fc_as_control_list(dots[[nm]], nm)
    arguments[[nm]] <- utils::modifyList(old, new)
    dots[[nm]] <- NULL
  }

  merged <- utils::modifyList(arguments, dots)
  keep <- names(merged) %in% names(formals(fun))
  merged[keep]
}


# ===========================================================================
# update.MIRT
# ===========================================================================

#' @describeIn update Update method for \code{MIRT} objects.
#'   Re-fits the MIRT model with modified arguments. Supports changing
#'   \code{model} (e.g., m2pl -> m3pl), \code{method} (iStEM -> Stan),
#'   \code{D}, \code{Q.matrix}, \code{control.model}, and
#'   \code{control.method}.
#' @method update MIRT
#' @export
update.MIRT <- function(object, ...) {
  args <- .fc_update_args(object$arguments, fit.MIRT, ...)
  do.call(fit.MIRT, args)
}

# ===========================================================================
# update.MGPCM
# ===========================================================================

#' @describeIn update Update method for \code{MGPCM} objects.
#'   Re-fits the MGPCM model with modified arguments.
#' @method update MGPCM
#' @export
update.MGPCM <- function(object, ...) {
  args <- .fc_update_args(object$arguments, fit.MGPCM, ...)
  do.call(fit.MGPCM, args)
}

# ===========================================================================
# update.MGGUM
# ===========================================================================

#' @describeIn update Update method for \code{MGGUM} objects.
#'   Re-fits the MGGUM model with modified arguments.
#' @method update MGGUM
#' @export
update.MGGUM <- function(object, ...) {
  args <- .fc_update_args(object$arguments, fit.MGGUM, ...)
  do.call(fit.MGGUM, args)
}

# ===========================================================================
# update.FCMIRT
# ===========================================================================

#' @describeIn update Update method for \code{FCMIRT} objects.
#'   Re-fits the FCMIRT model with modified arguments. Supports changing
#'   \code{model}, \code{method}, \code{fc.type}, \code{block.items},
#'   and control settings.
#' @method update FCMIRT
#' @export
update.FCMIRT <- function(object, ...) {
  args <- .fc_update_args(object$arguments, fit.FCMIRT, ...)
  do.call(fit.FCMIRT, args)
}

# ===========================================================================
# update.FCDCM
# ===========================================================================

#' @describeIn update Update method for \code{FCDCM} objects.
#'   Re-fits the FCDCM model with modified arguments. Supports changing
#'   \code{dcm.type}, \code{method}, \code{block.items}, and control
#'   settings.
#' @method update FCDCM
#' @export
update.FCDCM <- function(object, ...) {
  args <- .fc_update_args(object$arguments, fit.FCDCM, ...)
  do.call(fit.FCDCM, args)
}

# ===========================================================================
# update.FCGDINA
# ===========================================================================

#' @describeIn update Update method for \code{FCGDINA} objects.
#'   Re-fits the FCGDINA model with modified arguments.
#' @method update FCGDINA
#' @export
update.FCGDINA <- function(object, ...) {
  args <- .fc_update_args(object$arguments, fit.FCGDINA, ...)
  do.call(fit.FCGDINA, args)
}

# ===========================================================================
# update.FCGGUM
# ===========================================================================

#' @describeIn update Update method for \code{FCGGUM} objects.
#'   Re-fits the FCGGUM model with modified arguments.
#' @method update FCGGUM
#' @export
update.FCGGUM <- function(object, ...) {
  args <- .fc_update_args(object$arguments, fit.FCGGUM, ...)
  do.call(fit.FCGGUM, args)
}

# ===========================================================================
# update.TIRT
# ===========================================================================

#' @describeIn update Update method for \code{TIRT} objects.
#'   Re-fits the TIRT model with modified arguments. Supports changing
#'   \code{fc.type}, \code{method}, \code{pairs.value}, and control
#'   settings.
#' @method update TIRT
#' @export
update.TIRT <- function(object, ...) {
  args <- .fc_update_args(object$arguments, fit.TIRT, ...)
  do.call(fit.TIRT, args)
}
