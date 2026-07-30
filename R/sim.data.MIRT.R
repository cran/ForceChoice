#' Simulate Data from the Multidimensional IRT (MIRT) Model
#'
#' @description
#' Generates binary response data, latent trait vectors, and item parameters
#' from the multidimensional extension of the 1PL--4PL item response models.
#' The latent traits are drawn from a multivariate normal distribution with
#' a user-specified correlation matrix. Discrimination parameters respect a
#' Q-matrix structure and may be optionally rotated.
#'
#' @section Data Generation Process:
#'
#' \enumerate{
#'   \item \strong{Latent traits:}
#'         \eqn{\boldsymbol{\theta}_j \sim N_D(\mathbf{0}, \boldsymbol{\Sigma})},
#'         where \eqn{\boldsymbol{\Sigma}} is the \eqn{D \times D} correlation
#'         matrix \code{Corr}.
#'   \item \strong{Item parameters:}
#'         \itemize{
#'           \item \eqn{a_{id} \sim \text{Lognormal}(0.25, 0.25)} for
#'                 \eqn{q_{id} = 1} (M2PL/M3PL/M4PL);
#'                 \eqn{a_{id} = 1} for all items and dimensions in M1PL.
#'           \item \eqn{b_i \sim N(0, 1)} (difficulty).
#'           \item \eqn{c_i \sim U(0, 0.35)} (lower asymptote; M3PL/M4PL only).
#'           \item \eqn{d_i \sim U(0.65, 1)} (upper asymptote; M4PL only).
#'         }
#'   \item \strong{Response probabilities:}
#'         \eqn{P_{ij} = P(Y_{ij} = 1 \mid \boldsymbol{\theta}_j)} per the
#'         specified MIRT model (see \code{\link{fit.MIRT}} for the IRF
#'         equations).
#'   \item \strong{Binary responses:}
#'         \eqn{Y_{ij} \sim \text{Bernoulli}(P_{ij})}.
#' }
#'
#' Optional rotation of the loading matrix \code{a} and corresponding
#' transformation of \code{theta} and \code{Corr} is applied before
#' response generation.
#'
#' @param N Integer; number of examinees (default: 500).
#' @param I Integer; number of items (default: 20).
#' @param D Integer; number of latent dimensions (default: 2).
#' @param model Character; one of \code{"m1pl"}, \code{"m2pl"} (default),
#'   \code{"m3pl"}, \code{"m4pl"}.
#' @param Q.matrix An optional \eqn{I \times D} 0/1 matrix. If \code{NULL},
#'   a default structure with triangular identification pattern is used.
#' @param Corr An optional \eqn{D \times D} correlation matrix. If
#'   \code{NULL}, the identity matrix is used (orthogonal factors).
#' @param rotate Optional rotation method name (passed to
#'   \pkg{GPArotation} or \code{stats::promax}). Applicable only when
#'   \eqn{D > 1} and all Q-matrix entries are non-zero.
#' @param promax_m Integer; power parameter for \code{"promax"} rotation
#'   (default: 4).
#'
#' @return An object of class \code{"data.MIRT"}, a list containing:
#' \describe{
#'   \item{\code{data}, \code{response}}{\eqn{N \times I} binary response
#'         matrix (identical).}
#'   \item{\code{theta}}{\eqn{N \times D} matrix of true latent trait
#'         values.}
#'   \item{\code{par}}{\eqn{I \times (D+3)} matrix of true item parameters
#'         (columns \code{a1..aD, b, c, d}).}
#'   \item{\code{probability}}{\eqn{N \times I} matrix of response
#'         probabilities.}
#'   \item{\code{Q.matrix}}{\eqn{I \times D} Q-matrix used.}
#'   \item{\code{Corr}}{\eqn{D \times D} correlation matrix.}
#'   \item{\code{model}, \code{N}, \code{I}, \code{D}}{Data-generating
#'         parameters.}
#' }
#'
#' @seealso
#' \code{\link{fit.MIRT}} for fitting the MIRT model,
#' \code{\link{model.MIRT}} for computing response probabilities.
#'
#' @examples
#' # Basic 2D 2PL simulation
#' sim <- sim.data.MIRT(N = 20, I = 6, D = 2, model = "m2pl")
#'
#' # With correlated factors and rotation
#' Corr <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' sim_rot <- sim.data.MIRT(N = 20, I = 6, D = 2, model = "m2pl",
#'                          Q.matrix = matrix(1, 6, 2),
#'                          Corr = Corr, rotate = "oblimin")
#'
#' # Fit the simulated data
#' fit <- fit.MIRT(
#'   sim$response, model = "m2pl", D = 2, method = "iStEM",
#'   control.method = list(
#'     vis = FALSE, seed = 123,
#'     M = 2, B = 2, burnin.maxitr = 2,
#'     maxitr = 3, eps1 = 10, eps2 = 10,
#'     estimate.se = FALSE)
#' )
#' cor(fit$theta$est, sim$theta)  # trait recovery
#'
#' @export
sim.data.MIRT <- function(N = 500, I = 20, D = 2, model = "m2pl",
                          Q.matrix = NULL, Corr = NULL, rotate = NULL,
                          promax_m = 4) {

  call <- match.call()

  N <- check_integer_scalar(N, "N")
  I <- check_integer_scalar(I, "I")
  D <- check_integer_scalar(D, "D")
  model <- resolve_model_type(model)
  if (model == "m1pl" && D > 1L) {
    stop("M1PL simulation is supported only for D = 1.", call. = FALSE)
  }

  Corr <- validate_corr_matrix(Corr, D)

  theta <- generate_theta_mvn(N, D, Corr)

  Q.matrix <- istem_prepare_01_q(
    I, D, Q.matrix, triangular = TRUE, require.row = TRUE
  )

  if (model %in% c("m2pl", "m3pl", "m4pl")) {
    a <- matrix(rlnorm(I * D, meanlog = 0.25, sdlog = 0.25), I, D) * Q.matrix
  } else {
    a <- matrix(1, I, D)
  }

  b <- rnorm(I, 0, 1)
  c <- runif(I, 0.00, 0.35)
  d <- runif(I, 0.65, 1.00)

  if (model == "m1pl") {
    a <- matrix(1, I, D)
    c <- rep(0.0, I)
    d <- rep(1.0, I)
  } else if (model == "m2pl") {
    c <- rep(0.0, I)
    d <- rep(1.0, I)
  } else if (model == "m3pl") {
    d <- rep(1.0, I)
  }

  if (!is.null(rotate) && D > 1 && model != "m1pl") {
    if (any(Q.matrix == 0)) {
      stop("Rotation is not allowed when Q.matrix contains zeros (confirmatory MIRT). ",
           "The Q.matrix specifies a confirmatory structure that rotation would destroy.", call. = FALSE)
    }

    rot <- apply_rotation_promax(a, theta, Corr, rotate, D, promax_m)
    a     <- rot$a
    theta <- rot$theta
    Corr  <- rot$Corr
  }

  if (D > 1) {
    colnames(a) <- paste0("a", seq_len(D))
  }

  par <- cbind(a, b, c, d)

  probability <- model.MIRT(theta = theta, par = par)
  response <- +(probability > matrix(runif(N * I), nrow = N, ncol = I))
  colnames(response) <- seq_len(I)

  data.obj <- list(
    data        = response,
    response    = response,
    theta       = theta,
    par         = par,
    Q.matrix    = Q.matrix,
    block.items = NULL,
    Corr        = Corr,
    model       = model,
    N           = N,
    I           = I,
    D           = D,
    probability = probability,
    call        = call,
    arguments   = list(
      N         = N,
      I         = I,
      D         = D,
      model     = model,
      Q.matrix  = Q.matrix,
      Corr      = Corr,
      rotate    = rotate,
      promax_m  = promax_m
    )
  )

  class(data.obj) <- "data.MIRT"

  return(data.obj)
}


#' Compute MIRT Item Response Probabilities
#'
#' @description
#' Computes the probability of a correct response for each person--item
#' combination under the multidimensional 1PL--4PL models. This is the
#' core IRF (item response function) used throughout the package.
#'
#' @param theta An \eqn{N \times D} matrix of latent trait values.
#' @param par An \eqn{I \times (D+3)} parameter matrix with columns
#'   \code{a1..aD, b, c, d}, where \code{b} is the item difficulty.
#'   The linear predictor is \eqn{\mathbf{a}_i'\boldsymbol{\theta}_j - b_i}.
#'
#' @return An \eqn{N \times I} matrix of response probabilities
#'   \eqn{P(Y_{ij} = 1 \mid \boldsymbol{\theta}_j)}.
#'
#' @seealso \code{\link{fit.MIRT}} for the full model specification.
#' @keywords internal
model.MIRT <- function(theta, par) {
  cpp_model_MIRT(theta, par)
}
