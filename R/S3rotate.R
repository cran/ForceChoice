#' @title S3 Methods: rotate
#'
#' @description
#' Applies factor rotation to the posterior mean parameter estimates of
#' fitted exploratory multidimensional models from the \pkg{ForceChoice}
#' package. Rotation reparameterises the latent space to improve
#' interpretability of factor loadings while preserving the model-implied
#' covariance structure. Standard errors and other posterior uncertainty
#' summaries are carried over from the original fit object without
#' recomputation. No per-iteration rotation or cross-iteration alignment is
#' performed.
#'
#' For GPArotation methods, theta is re-standardized to a standard normal
#' scale after rotation. For Promax rotation, no post-rotation
#' standardization is applied. Loadings and the factor correlation matrix
#' are updated accordingly to preserve model-implied covariances.
#'
#' Rotation is not permitted when the fitted \code{Q.matrix} contains zeros
#' because that confirmatory structure would not generally be preserved
#' under rotation.
#'
#' @param object A fitted model object of class \code{"MIRT"} (returned by
#'   \code{\link[ForceChoice]{fit.MIRT}}) or \code{"MGPCM"} (returned by
#'   \code{\link[ForceChoice]{fit.MGPCM}}).
#' @param method Character string specifying the rotation method. Supported
#'   methods include \code{"promax"}, \code{"Varimax"}, \code{"quartimax"},
#'   \code{"oblimin"}, \code{"quartimin"}, \code{"geominT"}, \code{"geominQ"},
#'   \code{"targetT"}, \code{"targetQ"}, \code{"pstT"}, \code{"pstQ"}, and
#'   others provided by the \pkg{GPArotation} package.
#' @param vis Logical; if \code{TRUE} (default), displays informational
#'   messages during rotation.
#' @param target.a Numeric matrix of target loadings for target rotations
#'   (\code{"targetT"}, \code{"targetQ"}, \code{"pstT"}, \code{"pstQ"}).
#' @param target.theta Numeric matrix of target latent traits for target
#'   rotations (\code{"targetT"}, \code{"targetQ"} only).
#' @param pst.W Numeric weight matrix for PST rotations
#'   (\code{"pstT"}, \code{"pstQ"}).
#' @param lp.p Numeric; p parameter for Lp rotation (default = 1).
#' @param lp.gpaiter Integer; max GPA iterations for Lp rotation
#'   (default = 5).
#' @param promax_m Numeric; power parameter for Promax rotation
#'   (default = 4).
#' @param ... Additional arguments passed to the \pkg{GPArotation} rotation
#'   function.
#'
#' @return
#' \describe{
#'   \item{\code{rotate.MIRT}}{An object of class \code{"RotateMIRT"}
#'     containing the rotated \code{"MIRT"} fit object (\code{object}),
#'     the matched call (\code{call}), and the rotation arguments
#'     (\code{arguments}). Standard errors in the fitted object are
#'     preserved without recomputation.}
#'   \item{\code{rotate.MGPCM}}{An object of class \code{"RotateMGPCM"}
#'     containing the rotated \code{"MGPCM"} fit object (\code{object}),
#'     the matched call (\code{call}), and the rotation arguments
#'     (\code{arguments}). Category intercepts are retained from the
#'     original fit.}
#' }
#'
#' @details
#' Only posterior means are rotated. For \code{MGPCM} objects, the first
#' \eqn{D} columns of \code{object$par$est} are treated as discrimination
#' parameters; all category intercepts are unchanged.
#'
#' @seealso \code{\link[ForceChoice]{fit.MIRT}},
#'   \code{\link[ForceChoice]{fit.MGPCM}}
#' @name rotate
NULL

#' @export
rotate <- function(object, ...) {
  UseMethod("rotate")
}
