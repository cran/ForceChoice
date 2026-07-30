#' Fit the Multidimensional Generalized Graded Unfolding Model (MGGUM)
#'
#' @description
#' Fits a multidimensional generalization of the generalized graded unfolding
#' model (GGUM) to polytomous response data. Unlike dominance (cumulative)
#' IRT models, the MGGUM belongs to the ideal-point (unfolding) model family
#' where the probability of endorsement is highest when the person and item
#' locations match, decreasing as the person moves away in either direction.
#' Two estimation backends are provided: Bayesian MCMC (Stan) and iStEM.
#'
#' @section Model Specification:
#'
#' Let \eqn{Y_{ij} \in \{0, 1, \dots, K_i - 1\}} be the categorical response
#' of person \eqn{j} to item \eqn{i}, with \eqn{K_i \ge 2} categories.
#' Let \eqn{\boldsymbol{\theta}_j} denote the \eqn{D}-dimensional latent
#' trait vector.
#'
#' The category response probability implemented by \code{model.MGGUM},
#' iStEM, and Stan is the distance-based MGGUM form. Define
#' \deqn{
#'   r_{ij} =
#'   \left[
#'     \sum_{d=1}^{D} a_{id}^2
#'     \left(\theta_{jd} - \delta_{id}\right)^2
#'   \right]^{1/2},
#'   \qquad
#'   S_i = \sum_{d=1}^{D} a_{id}.
#' }
#' Let \eqn{\tau_{i0}=0} and
#' \eqn{\psi_{ik}=S_i\sum_{v=0}^{k}\tau_{iv}}. With
#' \eqn{M_i=2K_i-1}, the probability of category
#' \eqn{k = 0,\dots,K_i-1} is
#' \deqn{
#'   P(Y_{ij} = k \mid \boldsymbol{\theta}_j) =
#'   \frac{
#'     \exp\{k r_{ij}-\psi_{ik}\} +
#'     \exp\{(M_i-k)r_{ij}-\psi_{ik}\}
#'   }{
#'     \sum_{r=0}^{K_i-1}
#'     \bigl[
#'       \exp\{r r_{ij}-\psi_{ir}\} +
#'       \exp\{(M_i-r)r_{ij}-\psi_{ir}\}
#'     \bigr]
#'   }.
#' }
#' Here \eqn{a_{id} \ge 0} are discrimination parameters,
#' \eqn{\delta_{id}} are item-location parameters (signed positive when
#' \eqn{q_{id} = 1},
#' negative when \eqn{q_{id} = -1}), and threshold parameters \eqn{\tau_{ik}}
#' satisfying \eqn{\tau_{i0} = 0} and
#' \eqn{\tau_{i1} < \tau_{i2} < \dots < \tau_{i,K_i-1} < 0}.
#'
#' \strong{Q-matrix sign convention:}
#' \describe{
#'   \item{\eqn{q_{id} = 1}}{Active dimension; \eqn{\delta_{id} > 0}
#'         (item located on the positive side of dimension \eqn{d}).}
#'   \item{\eqn{q_{id} = -1}}{Active dimension; \eqn{\delta_{id} < 0}
#'         (item located on the negative side of dimension \eqn{d}).}
#'   \item{\eqn{q_{id} = 0}}{Inactive dimension;
#'         \eqn{a_{id} = \delta_{id} = 0}.}
#' }
#'
#' The sign of \eqn{q_{id}} controls the item-location side, not the sign of
#' the slope. Both \eqn{a_{id}} and \eqn{|\delta_{id}|} reflect the item's
#' relevance to dimension \eqn{d}.
#'
#' @section Estimation Methods:
#'
#' \describe{
#'   \item{\strong{Stan} (\code{method = "stan"}):}{
#'     Full Bayesian inference via HMC. All parameters are jointly sampled
#'     from the posterior distribution.}
#'   \item{\strong{iStEM} (\code{method = "iStEM"}):}{
#'     Improved Stochastic EM with finite-grid block Gibbs person-sampling
#'     and item optimization via L-BFGS-B.
#'     The threshold parameters \eqn{\tau_{ik}} are constrained to maintain
#'     monotonicity during optimization.}
#' }
#'
#' @param data An \eqn{N \times I} matrix of integer responses
#'   coded \eqn{0, 1, \dots, K_i-1}.
#' @param D Integer; number of latent dimensions (\eqn{D \ge 1}).
#'   Default is \code{2}.
#' @param Q.matrix An optional \eqn{I \times D} matrix with entries
#'   \eqn{-1, 0, 1}. See section \strong{Q-matrix sign convention}.
#'   Default generates a random signed matrix.
#' @param length.poly Integer vector (length \eqn{I}) or scalar giving the
#'   number of categories per item. If \code{NULL}, inferred from the
#'   observed response maxima.
#' @param method Estimation method: \code{"iStEM"} (default) or
#'   \code{"stan"}.
#' @param control.model A named list of model-level hyperparameters for the
#'   GGUM unfolding model. Supported entries:
#'   \describe{
#'     \item{\code{a.mu}, \code{a.sigma}}{Prior location and scale for
#'           \eqn{\log a_{id}} (log-normal). Defaults: 0, 0.5.
#'           The log-normal prior ensures positivity of discrimination
#'           parameters. The zero mean-log encourages slopes near 1
#'           unless the data strongly indicate otherwise.}
#'     \item{\code{delta.pos.mu}, \code{delta.pos.sigma}}{Prior mean and
#'           SD for \eqn{\delta_{id}} when \eqn{q_{id} = 1} (positive-side
#'           item locations; normal). Defaults: 1, 0.5. These priors
#'           regularize item locations on the positive side of each
#'           dimension.}
#'     \item{\code{delta.neg.mu}, \code{delta.neg.sigma}}{Prior mean and
#'           SD for \eqn{\delta_{id}} when \eqn{q_{id} = -1} (negative-side
#'           item locations; normal). Defaults: -1, 0.5. These priors
#'           regularize item locations on the negative side of each
#'           dimension.}
#'     \item{\code{tau.mu}, \code{tau.sigma}}{Location and scale for the
#'           log-normal prior on positive threshold magnitudes
#'           \eqn{-\tau_{ik}} for \eqn{k \ge 1}. Defaults: 0, 0.5.
#'           The thresholds satisfy
#'           \eqn{\tau_{i1} < \tau_{i2} < \dots < \tau_{i,K_i-1} < 0}
#'           with \eqn{\tau_{i0} = 0} fixed for identification.}
#'     \item{\code{theta.mu}}{Prior mean vector for
#'           \eqn{\boldsymbol{\theta}_j}. Default: \code{rep(0, D)}.}
#'     \item{\code{L}}{Theta grid size per dimension for marginal
#'           log-likelihood computation and iStEM block Gibbs sampling.
#'           Default adapts to \eqn{D}.}
#'     \item{\code{theta.lower}, \code{theta.upper}}{Bounds for the
#'           theta grid used by marginal log-likelihood computation and
#'           iStEM block Gibbs sampling. Defaults: -6, 6.}
#'   }
#' @param control.method A named list of method-specific tuning parameters.
#'   Common entries (used by both Stan and iStEM):
#'   \describe{
#'     \item{\code{cores}}{Number of CPU cores for parallel chains
#'           (Stan) or ignored (iStEM). Default: the number of
#'           \code{chains}.}
#'     \item{\code{vis}}{Logical; if \code{TRUE} (default), prints progress
#'           information to the console.}
#'     \item{\code{seed}}{Random seed for reproducibility.
#'           Default: a random integer.}
#'   }
#'   Stan-specific entries:
#'   \describe{
#'     \item{\code{chains}}{Number of MCMC chains (default: 2).}
#'     \item{\code{iter}}{Total iterations per chain (default: 5000).}
#'     \item{\code{warmup}}{Warmup/burn-in iterations per chain
#'           (default: \code{iter / 2}).}
#'     \item{\code{thin}}{Thinning interval (default: 1).}
#'     \item{\code{init}}{Initial values: \code{"random"} (default)
#'           for uniform(-2, 2) initialization, or a list of initial
#'           values per chain.}
#'     \item{\code{algorithm}}{MCMC algorithm: \code{"HMC"} (default),
#'           \code{"HMC"}, or \code{"Fixed_param"}.}
#'     \item{\code{adapt_delta}}{Target average acceptance probability
#'           (NUTS; default: 0.95). Values closer to 1 reduce step size
#'           and improve sampling for difficult posteriors.}
#'     \item{\code{max_treedepth}}{Maximum tree depth (NUTS; default: 10).
#'           Increase if "max treedepth exceeded" warnings appear.}
#'     \item{\code{stepsize}}{Initial step size for the leapfrog
#'           integrator (auto-tuned by Stan if not set).}
#'     \item{\code{int_time}}{Total integration time for HMC trajectories
#'           (only when \code{algorithm = "HMC"}).}
#'     \item{\code{metric}}{Mass matrix type: \code{"unit_e"},
#'           \code{"diag_e"} (default), or \code{"dense_e"}.}
#'     \item{\code{adapt_engaged}}{Logical; if \code{TRUE} (default),
#'           warmup adaptation is enabled.}
#'     \item{\code{adapt_init_buffer}, \code{adapt_term_buffer},
#'           \code{adapt_window}}{Warmup adaptation scheduling parameters
#'           (defaults: 25, 50, 25).}
#'   }
#'   iStEM-specific entries:
#'   \describe{
#'     \item{\code{M}}{Number of burn-in batches retained for Geweke
#'           convergence diagnosis (default: 10; must be \eqn{\ge 2}).}
#'     \item{\code{B}}{Batch size: MCMC iterations per batch (default: 20).}
#'     \item{\code{burnin.maxitr}}{Maximum burn-in batches (default: 100).}
#'     \item{\code{maxitr}}{Maximum total batches (default: 2000).}
#'     \item{\code{eps1}}{Geweke z-score convergence threshold
#'           (default: 1.5). The burn-in phase ends when
#'           \eqn{\sum z^2 / K < \epsilon_1} or the MC error criterion
#'           is satisfied.}
#'     \item{\code{eps2}}{Monte Carlo error tolerance (default: 0.4).
#'           The algorithm continues until
#'           \eqn{max(d_k \cdot N) < \epsilon_2}.}
#'     \item{\code{frac1}, \code{frac2}}{Fractions for the Geweke
#'           diagnostic (defaults: 0.1, 0.5).}
#'     \item{\code{corr.optim.maxit}}{Maximum L-BFGS-B iterations for the
#'           constrained unit-diagonal correlation update (default: 50).}
#'     \item{\code{optim.maxit}}{Maximum L-BFGS-B iterations per item
#'           during the item-parameter update (default: 50). Increase
#'           if item optimization does not converge within this limit.}
#'     \item{\code{fix.corr}}{Logical; if \code{TRUE}, the inter-trait
#'           correlation matrix is fixed to the identity (default:
#'           \code{FALSE}).}
#'     \item{\code{estimate.se}}{Logical; if \code{TRUE} (default),
#'           standard errors are computed from the final Monte Carlo
#'           chain.}
#'     \item{\code{delta.lower}}{Lower bound on \eqn{|\delta_{id}|}
#'           (default: \code{1e-4}). Ensures item locations are bounded
#'           away from zero for active dimensions.}
#'     \item{\code{tau.gap.lower}, \code{tau.gap.upper}}{Lower and upper
#'           bounds on the positive gaps used to construct ordered
#'           thresholds. Defaults are 1e-4 and 6.}
#'   }
#'
#' @return An object of class \code{"MGGUM"} containing:
#' \describe{
#'   \item{\code{npar}}{Number of free parameters.}
#'   \item{\code{method}}{\code{"stan"} or \code{"iStEM"}.}
#'   \item{\code{theta}}{List with \code{est}, \code{se}, \code{Rhat}
#'         (\eqn{N \times D}).}
#'   \item{\code{par}}{List with \code{est}, \code{se}, \code{Rhat},
#'         \code{free} (\eqn{I \times (2D + K_{max})}).
#'         Columns: \code{a1..aD, delta1..deltaD, tau0..tau_{K_{max}-1}}.}
#'   \item{\code{Corr}}{List with \code{est}, \code{se}, \code{Rhat}
#'         (\eqn{D \times D}).}
#'   \item{\code{length.poly}}{Per-item category counts.}
#'   \item{\code{logLik}}{Marginal log-likelihood (class \code{"logLik"}).}
#' }
#'
#' @references
#' Roberts, J. S., Donoghue, J. R., & Laughlin, J. E. (2000). A general item
#'   response theory model for unfolding unidimensional polytomous responses.
#'   \emph{Applied Psychological Measurement}, 24(1), 3--32.
#'   \doi{10.1177/01466216000241001}
#'
#' Usami, S. (2011). Generalized graded unfolding model with structural
#'   equation for subject parameters. \emph{Japanese Psychological Research},
#'   53(3), 221--232.
#'
#' @seealso
#' \code{\link{sim.data.MGGUM}}, \code{\link{get.fit.index.MGGUM}},
#' \code{\link{logLik.MGGUM}}
#'
#' @examples
#' sim <- sim.data.MGGUM(N = 20, I = 6, D = 2, length.poly = 4)
#' fit <- fit.MGGUM(sim$response, D = 2, method = "iStEM",
#'                  control.method = list(
#'                    vis = FALSE, seed = 123,
#'                    M = 2, B = 2, burnin.maxitr = 2,
#'                    maxitr = 3, eps1 = 10, eps2 = 10,
#'                    estimate.se = FALSE))
#' head(fit$theta$est)
#' gof <- get.fit.index(fit)
#' summary(gof)
#'
#' @export
fit.MGGUM <- function(data, D = NULL, Q.matrix = NULL, length.poly = NULL,
                      method = c("iStEM", "stan"),
                      control.model = NULL,
                      control.method = NULL) {

  call <- match.call()
  method <- match.arg(method)
  data <- istem_prepare_response(data, binary = FALSE)
  length.poly <- istem_prepare_length_poly(data, length.poly)
  D <- resolve_D(D, Q.matrix = Q.matrix, I = ncol(data))
  Q.matrix <- istem_prepare_signed_q(
    I = ncol(data), D = D, Q.matrix = Q.matrix
  )
  control.model  <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")

  if (method == "stan") {
    return(fit.MGGUM.stan(
      response = data,
      D = D,
      Q.matrix = Q.matrix,
      length.poly = length.poly,
      control.model = control.model,
      control.method = control.method,
      .call = call
    ))
  }

  fit.MGGUM.iStEM(
    response = data,
    D = D,
    Q.matrix = Q.matrix,
    length.poly = length.poly,
    control.model = control.model,
    control.method = control.method,
    .call = call
  )
}
