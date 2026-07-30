#' Simulate Data from the Multidimensional Generalized Graded Unfolding Model
#'
#' @description
#' Generates polytomous responses from the same distance-based MGGUM
#' probability function used by \code{\link{model.MGGUM}} and
#' \code{\link{fit.MGGUM}}. For each item, active discriminations are drawn
#' from \eqn{U(0.5, 2)}, active locations follow the sign of the Q-matrix,
#' \eqn{\tau_{i0}=0}, and the remaining \eqn{\tau} values are ordered
#' negative thresholds.
#'
#' @section Data Generation Process:
#' \enumerate{
#'   \item Draw latent traits from
#'         \eqn{N_D(\mathbf{0}, \boldsymbol{\Sigma})}.
#'   \item Generate discrimination, location, and threshold parameters under
#'         the sign constraints encoded by \code{Q.matrix}.
#'   \item Compute stacked category probabilities with
#'         \code{\link{model.MGGUM}}.
#'   \item Sample one ordinal response per person and item. A category can be
#'         absent in a finite sample; \code{length.poly} retains the intended
#'         support.
#' }
#'
#' @param N Integer; number of examinees (default: 500).
#' @param I Integer; number of items (default: 20).
#' @param D Integer; number of latent dimensions (default: 2).
#' @param length.poly Integer vector or scalar indicating the number of
#'   categories for each item. A scalar is recycled to length \eqn{I}.
#' @param Q.matrix Optional \eqn{I \times D} matrix with values -1, 0, or 1.
#'   \itemize{
#'     \item \code{1}: active dimension with positive-side \code{delta}
#'     \item \code{-1}: active dimension with negative-side \code{delta}
#'     \item \code{0}: \code{a} and \code{delta} fixed to 0 (inactive dimension)
#'   }
#'   The sign controls the item-location side, not the sign of discrimination.
#'   If \code{NULL}, a default matrix of randomly signed active entries is used.
#' @param Corr Optional \eqn{D \times D} correlation matrix. If \code{NULL},
#'   the identity matrix is used.
#'
#' @return An object of class \code{"data.MGGUM"}, a list containing:
#' \describe{
#'   \item{\code{data}, \code{response}}{\eqn{N \times I} integer response
#'         matrices with categories coded from 0 to \eqn{K_i - 1}.}
#'   \item{\code{theta}}{\eqn{N \times D} matrix of true latent traits.}
#'   \item{\code{par}}{\eqn{I \times (2D + \max_i K_i)} matrix of GGUM item
#'         parameters: discrimination columns, location columns, and stacked
#'         threshold columns.}
#'   \item{\code{probability}}{\eqn{N \times \sum_i K_i} matrix of stacked
#'         category probabilities.}
#'   \item{\code{Q.matrix}, \code{length.poly}, \code{Corr}}{Design matrix,
#'         category counts, and latent correlation matrix used to generate
#'         the data.}
#'   \item{\code{N}, \code{I}, \code{D}, \code{call}}{Data-generating
#'         metadata.}
#' }
#'
#' @seealso \code{\link{fit.MGGUM}}, \code{\link{model.MGGUM}}
#'
#' @examples
#' set.seed(123)
#' sim <- sim.data.MGGUM(N = 20, I = 5, D = 2, length.poly = 4)
#' str(sim$response)
#' dim(sim$par)
#' table(sim$response[, 1])
#'
#' @export
sim.data.MGGUM <- function(N = 500, I = 20, D = 2, length.poly = 5,
                           Q.matrix = NULL, Corr = NULL) {
  call <- match.call()

  N <- check_integer_scalar(N, "N")
  I <- check_integer_scalar(I, "I")
  D <- check_integer_scalar(D, "D")
  length.poly <- normalize_length_poly(length.poly, I)
  max_poly <- max(length.poly)

  Corr <- validate_corr_matrix(Corr, D)

  if (is.null(Q.matrix)) {
    Q.matrix <- matrix(sample(c(-1, 1), I * D, replace = TRUE), I, D)
  }
  Q.matrix <- istem_prepare_signed_q(I, D, Q.matrix)

    theta <- generate_theta_mvn(N, D, Corr)

    # Generate a parameters respecting Q.matrix
    # a ~ U(0.50, 2.00), free when |Q| = 1, fixed to 0 when Q = 0
    a <- matrix(runif(I * D, 0.5, 2.0), I, D) * abs(Q.matrix)

    # par matrix: I x (D + D + max_poly)
    # Columns: a1..aD | delta1..deltaD | tau1(=0) tau2..tau_{max_poly}
    par <- matrix(NA, nrow = I, ncol = D + D + max_poly)
    par[, 1:D] <- a
    mask <- matrix(NA, nrow = I, ncol = D + D + max_poly)
    mask[, 1:D] <- 1
    for (i in seq_len(I)) {

      # delta parameters: sign determined by Q.matrix, magnitude ~ U(0, 2.00)
      for (d in seq_len(D)) {
        if (Q.matrix[i, d] == 1) {
          par[i, D + d] <- runif(1, 0.0, 2.0)
        } else if (Q.matrix[i, d] == -1) {
          par[i, D + d] <- runif(1, -2.0, 0.0)
        } else {
          par[i, D + d] <- 0.0
        }
      }
      mask[i, (D + 1):(D + D)] <- 1

      # tau parameters: tau0 = 0 (fixed), tau1 < tau2 < ... < tau_{K-1} < 0
      par[i, D + D + 1] <- 0.0
      mask[i, D + D + 1] <- 1
      Ki <- length.poly[i]
      if (Ki > 1) {
        Kf <- Ki - 1  # number of free tau parameters
        tau_free <- numeric(Kf)
        for (t in seq_len(Kf)) {
          lo <- -2.00 + (t - 1) * 2.0 / Kf
          hi <- -2.00 +  t      * 2.0 / Kf
          tau_free[t] <- runif(1, lo, hi)
        }
        # Sort ascending to ensure tau1 < tau2 < ... < tau_{K-1} < 0
        tau_sorted <- sort(tau_free, decreasing = FALSE)
        par[i, (D + D + 2):(D + D + Ki)] <- tau_sorted
        mask[i, (D + D + 2):(D + D + Ki)] <- 1
      }
    }

    prob <- model.MGGUM(theta = theta, par = par)
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
  colnames(par) <- c(paste0("a", seq_len(D)),
                     paste0("delta", seq_len(D)),
                     paste0("tau", 0:(max_poly - 1)))
  rownames(par) <- paste0("item", seq_len(I))
  colnames(response) <- paste0("item", seq_len(I))

  ak <- matrix(0:(max_poly - 1), I, max_poly, byrow = TRUE)
  colnames(ak) <- paste0("ak", 0:(max_poly - 1))
  par.mirt <- cbind(par[, 1:D, drop = FALSE],
                    par[, (D + 1):(D + D), drop = FALSE],
                    ak,
                    par[, (D + D + 1):ncol(par), drop = FALSE])

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
      Corr        = Corr
    )
  )

  class(data.obj) <- "data.MGGUM"

  return(data.obj)
}

#' Compute MGGUM Category Response Probabilities
#'
#' @param theta \eqn{N \times D} matrix of latent trait values.
#' @param par \eqn{I \times (2D + \max K_i)} parameter matrix with columns
#'   \code{a1..aD, delta1..deltaD, tau0, tau1, ...}. Missing category
#'   columns for shorter items should be \code{NA}; \code{tau0} is fixed at 0.
#' @return An \eqn{N \times \sum_i K_i} matrix of stacked category
#'   probabilities. The implementation uses
#'   \eqn{r_{ij} = [\sum_d a_{id}^2(\theta_{jd}-\delta_{id})^2]^{1/2}}
#'   and
#'   \eqn{P(Y_{ij}=k) \propto
#'   \exp\{k r_{ij}-\psi_{ik}\} +
#'   \exp\{(2K_i-1-k)r_{ij}-\psi_{ik}\}}.
#' @keywords internal
model.MGGUM <- function(theta, par){
  cpp_model_MGGUM(theta, par)
}
