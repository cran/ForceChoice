#' Simulate Data from the Multidimensional Generalized Partial Credit Model
#'
#' @description
#' Generates polytomous (multi-category) response data from the MGPCM.
#' Category probabilities use the same softmax formulation as
#' \code{\link{model.MGPCM}}: \eqn{\log P(Y=k)} is proportional to
#' \eqn{k\,\mathbf{a}_i'\boldsymbol{\theta}_j + d_{ik}}.
#' The requested category counts are retained even when a finite sample does
#' not contain every possible category.
#'
#' @section Data Generation Process:
#' \enumerate{
#'   \item \strong{Latent traits:}
#'         \eqn{\boldsymbol{\theta}_j \sim N_D(\mathbf{0}, \boldsymbol{\Sigma})}.
#'   \item \strong{Item parameters:}
#'         \eqn{a_{id} \sim \text{Lognormal}(0.25, 0.25)} (active dimensions),
#'         \eqn{d_{i0}=0}, and
#'         \eqn{d_{i1},\dots,d_{i,K_i-1}} are sorted draws from
#'         \eqn{N(0,1)} used as category intercepts.
#'   \item \strong{Responses:}
#'         \eqn{Y_{ij} \sim \text{Categorical}(P(Y_{ij} = k \mid \boldsymbol{\theta}_j))}
#'         for \eqn{k = 0, \dots, K_i - 1}.
#' }
#'
#' @param N Integer; number of examinees (default: 500).
#' @param I Integer; number of items (default: 20).
#' @param D Integer; number of latent dimensions (default: 2).
#' @param length.poly Integer vector or scalar; number of categories per
#'   item. Recycled to length \eqn{I} if scalar (default: 5).
#' @param Q.matrix Optional \eqn{I \times D} 0/1 matrix. \code{NULL} uses
#'   a triangular identification pattern.
#' @param Corr Optional \eqn{D \times D} correlation matrix.
#' @param rotate Optional rotation method (\pkg{GPArotation} or
#'   \code{"promax"}).
#' @param promax_m Power parameter for Promax (default: 4).
#'
#' @return An object of class \code{"data.MGPCM"}, a list containing:
#' \describe{
#'   \item{\code{data}, \code{response}}{\eqn{N \times I} integer response
#'         matrices with categories coded from 0 to \eqn{K_i - 1}.}
#'   \item{\code{theta}}{\eqn{N \times D} matrix of true latent traits.}
#'   \item{\code{par}}{\eqn{I \times (D + \max_i K_i)} matrix of item
#'         parameters; discrimination columns are followed by category
#'         intercept columns.}
#'   \item{\code{probability}}{\eqn{N \times \sum_i K_i} matrix of stacked
#'         category probabilities.}
#'   \item{\code{Q.matrix}, \code{length.poly}, \code{Corr}}{Design matrix,
#'         category counts, and latent correlation matrix used to generate
#'         the data.}
#'   \item{\code{model}, \code{N}, \code{I}, \code{D}, \code{call}}{
#'         Data-generating metadata.}
#' }
#'
#' @seealso \code{\link{fit.MGPCM}}, \code{\link{model.MGPCM}}
#'
#' @examples
#' set.seed(123)
#' sim <- sim.data.MGPCM(N = 20, I = 5, D = 2, length.poly = 4)
#' str(sim$response)
#' dim(sim$probability)
#' sim$length.poly
#'
#' @export
sim.data.MGPCM <- function(N = 500, I = 20, D = 2, length.poly = 5,
                           Q.matrix = NULL, Corr = NULL, rotate = NULL,
                           promax_m = 4) {
  call <- match.call()

  N <- check_integer_scalar(N, "N")
  I <- check_integer_scalar(I, "I")
  D <- check_integer_scalar(D, "D")
  length.poly <- normalize_length_poly(length.poly, I)
  max_poly <- max(length.poly)

  Corr <- validate_corr_matrix(Corr, D)

  Q.matrix <- istem_prepare_01_q(I, D, Q.matrix, triangular = TRUE)

    theta <- generate_theta_mvn(N, D, Corr)

    a <- matrix(rlnorm(I * D, meanlog = 0.25, sdlog = 0.25), I, D) * Q.matrix

    if (!is.null(rotate) && D > 1) {
      rot <- apply_rotation_promax(a, theta, Corr, rotate, D, promax_m)
      a     <- rot$a
      theta <- rot$theta
      Corr  <- rot$Corr
    }

    par <- matrix(NA, nrow = I, ncol = D + max_poly)
    par[, 1:D] <- a
    mask <- matrix(NA, nrow = I, ncol = D + max_poly)
    mask[, 1:D] <- 1
    for (i in seq_len(I)) {

      par[i, D + 1] <- 0.0
      mask[i, D + 1] <- 1
      if (length.poly[i] > 1) {
        par[i, (D + 2):(D + length.poly[i])] <- sort(rnorm(length.poly[i] - 1, 0, 1))
        mask[i, (D + 2):(D + length.poly[i])] <- 1
      }
    }

    prob <- model.MGPCM(theta = theta, par = par)
    response <- matrix(0, N, I)
    idx <- 0
    for (i in seq_len(I)) {
      P.i <- prob[, (idx + 1):(idx + length.poly[i]), drop = FALSE]
      idx <- idx + length.poly[i]
      cdf <- t(apply(P.i, 1L, cumsum))
      response[, i] <- pmin(
        rowSums(cdf < stats::runif(N)), length.poly[i] - 1L
      )
    }

  if (D > 1) {
    colnames(a) <- paste0("a", seq_len(D))
  }
  colnames(par) <- c(paste0("a", seq_len(D)), paste0("d", 0:(max_poly - 1)))
  rownames(par) <- paste0("item", seq_len(I))
  colnames(response) <- paste0("item", seq_len(I))

  ak <- matrix(0:(max_poly - 1), I, max_poly, byrow = TRUE)
  colnames(ak) <- paste0("ak", 0:(max_poly - 1))
  par.mirt <- cbind(par[, 1:D, drop = FALSE], ak,
                    par[, (D + 1):ncol(par), drop = FALSE])

  data.obj <- list(
    data        = response,
    response    = response,
    theta       = theta,
    par         = par,
    par.mirt    = par.mirt,
    Q.matrix    = Q.matrix,
    block.items = NULL,
    Corr        = Corr,
    length.poly = length.poly,
    mask        = mask,
    probability = prob,
    N           = N,
    I           = I,
    D           = D,
    call        = call,
    arguments   = list(
      N           = N,
      I           = I,
      D           = D,
      length.poly = length.poly,
      Q.matrix    = Q.matrix,
      Corr        = Corr,
      rotate      = rotate,
      promax_m    = promax_m
    )
  )

  class(data.obj) <- "data.MGPCM"

  return(data.obj)
}

#' Compute MGPCM Category Response Probabilities
#'
#' @param theta \eqn{N \times D} matrix of latent trait values.
#' @param par \eqn{I \times (D + \max K_i)} parameter matrix with
#'   columns \code{a1..aD, d0, d1, ...}. Missing category columns for
#'   shorter items should be \code{NA}; \code{d0} is typically fixed at 0.
#' @return An \eqn{N \times \sum_i K_i} matrix of category probabilities
#'   (stacked by item), where
#'   \eqn{P(Y_{ij}=k) \propto \exp\{k\,\mathbf{a}_i'
#'   \boldsymbol{\theta}_j + d_{ik}\}}.
#' @keywords internal
model.MGPCM <- function(theta, par){
  cpp_model_MGPCM(theta, par)
}
