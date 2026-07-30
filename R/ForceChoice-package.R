#' ForceChoice: Forced-Choice Modeling Based on IRT and CDM
#'
#' @description
#' The \pkg{ForceChoice} package provides a unified framework for fitting,
#' simulating, and evaluating forced-choice and traditional item response
#' models. It supports eight model families spanning dominance (MIRT),
#' ideal-point/unfolding (MGGUM), polytomous (MGPCM), and forced-choice
#' (FCMIRT, FCGGUM, TIRT, FCDCM, FCGDINA) paradigms, with full Bayesian
#' estimation via Stan, fast improved stochastic EM (iStEM), and a
#' deterministic EM backend for FCGDINA.
#'
#' @section Model Families:
#'
#' \strong{Traditional (single-stimulus) models:}
#' \describe{
#'   \item{\code{\link{fit.MIRT}}}{Multidimensional IRT (1PL--4PL).
#'         Binary responses with optional Q-matrix structure.}
#'   \item{\code{\link{fit.MGPCM}}}{Multidimensional Generalized Partial
#'         Credit Model. Polytomous ordered-category responses.}
#'   \item{\code{\link{fit.MGGUM}}}{Multidimensional Generalized Graded
#'         Unfolding Model. Ideal-point polytomous responses.}
#' }
#'
#' \strong{Forced-choice (comparative) models:}
#' \describe{
#'   \item{\code{\link{fit.FCMIRT}}}{Forced-choice MIRT. Item-level
#'         dominance model with sequential ranking over item endorsement
#'         logits at the block level.}
#'   \item{\code{\link{fit.FCGGUM}}}{Forced-choice GGUM. Item-level
#'         ideal-point model with sequential ranking over binary GGUM
#'         endorsement logits.}
#'   \item{\code{\link{fit.TIRT}}}{Thurstonian IRT for forced-choice.
#'         Pairwise probit comparison of latent utility differences.}
#'   \item{\code{\link{fit.FCDCM}}}{Forced-choice DCM. Higher-order
#'         cognitive diagnostic model with exact attribute-profile
#'         marginalization.}
#'   \item{\code{\link{fit.FCGDINA}}}{Forced-choice GDINA. General CDM
#'         item-response structure with forced-choice block likelihood.}
#' }
#'
#' @section Estimation Methods:
#'
#' Model families support a common estimation interface:
#' \describe{
#'   \item{\strong{Stan} (\code{method = "stan"}):}{
#'         Full Bayesian inference via Hamiltonian Monte Carlo
#'         (NUTS/HMC). Provides posterior means, standard deviations,
#'         and convergence diagnostics (\eqn{\hat{R}}).}
#'   \item{\strong{iStEM} (\code{method = "iStEM"}):}{
#'         Improved Stochastic EM (iStEM) algorithm
#'         combining finite-grid block Gibbs person sampling with
#'         L-BFGS-B item-parameter optimization. Scales well to moderate
#'         data sets.}
#'   \item{\strong{EM} (\code{method = "EM"}):}{
#'         Deterministic posterior-weight EM for FCGDINA. Useful as a
#'         reproducible baseline or fast diagnostic estimator.}
#' }
#'
#' @section Core Workflow:
#'
#' A typical analysis follows four steps:
#' \enumerate{
#'   \item \strong{Simulate} (or load) data with \code{sim.data.*()}.
#'   \item \strong{Fit} a model with \code{fit.*()}.
#'   \item \strong{Rotate} the solution if needed with
#'         the \code{\link{rotate}} S3 generic.
#'   \item \strong{Evaluate} fit with \code{\link{get.fit.index}()}
#'         and \code{\link{summary.good.of.fit}()}.
#' }
#'
#' @section Goodness-of-Fit:
#'
#' The package uses the limited-information M2 framework
#' (Maydeu-Olivares & Joe, 2005, 2006) implemented in
#' \code{\link{good.of.fit}}. Fit indices include:
#' \itemize{
#'   \item Information criteria: AIC, AICc, BIC, CAIC, SABIC, HQIC
#'   \item Absolute fit: M2, RMSEA (with CI), SRMSR, McDonald NCI
#'   \item Comparative fit: CFI, TLI, IFI
#'   \item Pseudo-\eqn{R^2}: McFadden, Cox--Snell, Nagelkerke, etc.
#'   \item Local dependence: Yen's Q3 (raw and adjusted)
#'   \item Classification: posterior entropy, mean max posterior
#' }
#'
#' @section Reproducibility:
#'
#' The package is designed so that simulation, estimation, extraction, and
#' model checking can be scripted end to end. Simulation functions return the
#' true parameters and the effective call arguments; fitted objects retain the
#' original input arguments, estimates, convergence diagnostics, and
#' likelihood-related quantities used by \code{\link{logLik}} and
#' \code{\link{get.fit.index}}. For computationally expensive Stan or iStEM
#' workflows, manuscript replication code should set seeds explicitly and use
#' a small demonstration configuration in addition to any full-scale analysis.
#'
#' @examples
#' set.seed(123)
#' sim <- sim.data.FCGDINA(N.person = 12, N.block = 2, I.block = 2,
#'                         D = 2, model = "DINA")
#' str(sim$response)
#' head(sim$Q.matrix)
#'
#' @name ForceChoice-package
#' @aliases ForceChoice-package
#' @useDynLib ForceChoice, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom rstan sampling
#' @import rstantools
#' @importFrom stats fitted plogis pnorm qlogis qnorm rlnorm rnorm runif sd setNames
#' @importFrom utils combn flush.console tail
"_PACKAGE"
