#' @title Extract Log-Likelihood from Fitted Models
#'
#' @description
#' Extracts the marginal log-likelihood from a fitted ForceChoice model
#' object. For continuous-trait IRT models, the log-likelihood is computed
#' by numerically integrating the joint likelihood of the observed responses
#' over the latent trait distribution using Cartesian-product normal
#' quadrature. For FCGDINA, the marginal likelihood is obtained by summing
#' over the \eqn{2^D} discrete attribute profiles.
#'
#' For a sample of \eqn{N} independent persons with responses
#' \eqn{\mathbf{Y}_j}, the marginal log-likelihood is:
#' \deqn{
#'   \ell = \sum_{j=1}^{N} \log \int
#'   P(\mathbf{Y}_j \mid \boldsymbol{\theta})\,
#'   \phi_D(\boldsymbol{\theta}; \mathbf{0}, \boldsymbol{\Sigma})\,
#'   d\boldsymbol{\theta},
#' }
#' where \eqn{\phi_D(\cdot)} is the \eqn{D}-variate normal density with
#' correlation matrix \eqn{\boldsymbol{\Sigma}}. The integral is approximated
#' on a Cartesian-product quadrature grid over \eqn{[-6, 6]^D}.
#'
#' The continuous-trait quadrature weights are normalized normal densities
#' computed on the grid points, with the fitted correlation matrix used when
#' available. FCGDINA instead uses the estimated class-proportion vector
#' \eqn{\boldsymbol{\pi}} over attribute profiles.
#'
#' @param object A fitted model object of class \code{"MIRT"}, \code{"MGPCM"},
#'   \code{"MGGUM"}, \code{"FCMIRT"}, \code{"FCDCM"}, \code{"FCGDINA"},
#'   \code{"FCGGUM"}, or \code{"TIRT"}.
#' @param theta.low Numeric; lower bound for the quadrature grid on the
#'   latent scale (default: \code{-6}).
#' @param theta.up Numeric; upper bound for the quadrature grid (default:
#'   \code{6}).
#' @param L Integer; number of quadrature points per dimension. If
#'   \code{NULL}, automatically chosen based on \eqn{D} (e.g., 61 for
#'   \eqn{D = 1}, 31 for \eqn{D = 2}, 15 for \eqn{D = 3}).
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"logLik"} with attributes:
#' \itemize{
#'   \item \code{nobs}: number of observations (\eqn{N})
#'   \item \code{df}: number of free parameters
#'   \item \code{theta.norm}: continuous-trait quadrature grid points, when
#'         applicable
#'   \item \code{prob}: model-implied response probabilities at grid points,
#'         when applicable
#'   \item \code{pi}: quadrature weights or FCGDINA class probabilities
#'   \item \code{L}: quadrature grid size used, when applicable
#'   \item \code{prob.class}: FCGDINA class-conditional forced-choice
#'         probabilities, when applicable
#' }
#'
#' @seealso
#' \code{\link{fit.MIRT}}, \code{\link{get.fit.index}},
#' \code{\link{good.of.fit}}
#'
#' @examples
#' sim <- sim.data.MIRT(N = 20, I = 6, D = 2, model = "m2pl")
#' fit <- fit.MIRT(
#'   sim$response, model = "m2pl", D = 2, method = "iStEM",
#'   control.method = list(
#'     vis = FALSE, seed = 123,
#'     M = 2, B = 2, burnin.maxitr = 2,
#'     maxitr = 3, eps1 = 10, eps2 = 10,
#'     estimate.se = FALSE)
#' )
#' ll <- logLik(fit, L = 9)
#' as.numeric(ll)
#' attr(ll, "df")
#'
#' @name logLik-class
#' @rdname logLik-class
#' @importFrom stats logLik
NULL

# The following suppresses roxygen2 auto-detection of class<-"logLik"
# as an S4 class definition.
#' @title S3 Methods for logLik
#'
#' @description
#' S3 methods for extracting marginal log-likelihoods from fitted
#' \pkg{ForceChoice} model objects.
#'
#' @param object A fitted model object of class \code{"MIRT"},
#'   \code{"MGPCM"}, \code{"MGGUM"}, \code{"FCMIRT"}, \code{"FCDCM"},
#'   \code{"FCGDINA"}, \code{"FCGGUM"}, or \code{"TIRT"}.
#' @param theta.low Numeric; lower bound for the quadrature grid on the
#'   latent scale. Used by continuous-trait models; if \code{NULL}, the value
#'   stored in the fitted object's controls is used, falling back to
#'   \code{-6}.
#' @param theta.up Numeric; upper bound for the quadrature grid on the latent
#'   scale. Used by continuous-trait models; if \code{NULL}, the value stored
#'   in the fitted object's controls is used, falling back to \code{6}.
#' @param L Integer; number of quadrature points per dimension. If
#'   \code{NULL}, the package chooses a dimension-dependent default.
#' @param ... Additional arguments passed to the method. Currently ignored by
#'   the implemented methods.
#'
#' @return An object of class \code{"logLik"}: a length-one numeric vector
#'   containing the marginal log-likelihood, with at least attributes
#'   \code{nobs} (sample size) and \code{df} (number of free parameters).
#'   Continuous-trait methods also attach quadrature diagnostics such as
#'   \code{theta.norm}, \code{prob}, \code{pi}, and \code{L}; FCGDINA attaches
#'   latent-class probabilities and class-conditional response probabilities.
#'
#' @details
#' For MIRT, MGPCM, MGGUM, FCMIRT, FCGGUM, TIRT, and FCDCM objects, the
#' likelihood is evaluated by marginalizing over a normal quadrature grid
#' (and, for FCDCM, over attribute profiles conditional on the higher-order
#' trait). For FCGDINA objects, the likelihood is evaluated exactly over the
#' \eqn{2^D} latent attribute profiles.
#'
#' @name logLik-class
#' @rdname logLik-class
#' @keywords internal
NULL
