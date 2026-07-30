#' Fit the Forced-Choice Generalized Graded Unfolding Model (FCGGUM)
#'
#' @description
#' Fits a forced-choice variant of the multidimensional generalized graded
#' unfolding model (GGUM) to comparative ranking data. Unlike dominance-based
#' FC models, the FCGGUM uses an ideal-point (unfolding) item response
#' process: endorsement probability peaks when the person and item locations
#' coincide, and the block-level ranking follows the Luce--Plackett model.
#'
#' @section Model Specification:
#'
#' \strong{Item-level unfolding process.}
#' For statement \eqn{i} with Q-vector \eqn{\mathbf{q}_i \in \{-1, 0, 1\}^D},
#' the category response function is the GGUM (see \code{\link{fit.MGGUM}}
#' for the full polytomous specification). For the forced-choice context,
#' each item is binary (\eqn{K_i=2}). Define
#' \deqn{
#'   r_{ij} =
#'   \left[
#'     \sum_{d=1}^{D} a_{id}^2
#'     \left(\theta_{jd} - \delta_{id}\right)^2
#'   \right]^{1/2},
#'   \qquad
#'   S_i = \sum_{d=1}^{D} a_{id},
#'   \qquad
#'   \psi_{i1} = \tau_{i1}S_i .
#' }
#' The binary GGUM endorsement probability used by the code is:
#' \deqn{
#'   P_i(\boldsymbol{\theta}_j) =
#'   \frac{
#'     \exp\{r_{ij}-\psi_{i1}\} +
#'     \exp\{2r_{ij}-\psi_{i1}\}
#'   }{
#'     1 + \exp\{3r_{ij}\} +
#'     \exp\{r_{ij}-\psi_{i1}\} +
#'     \exp\{2r_{ij}-\psi_{i1}\}
#'   },
#' }
#' where \eqn{\tau_{i1} < 0} is the first (and only free) threshold.
#'
#' \strong{Block-level forced-choice ranking.}
#' Within each block \eqn{b}, the probability of the observed ranking is
#' computed by applying the sequential Luce/Plackett rule to
#' \eqn{u_i(\boldsymbol{\theta}_j)=\mathrm{logit}
#' \{P_i(\boldsymbol{\theta}_j)\}}:
#' \deqn{
#'   P(\text{ranking} \mid \boldsymbol{\theta}_j) =
#'   \prod_{m=1}^{K_b-1}
#'   \frac{\exp\{u_{i_{(m)}}(\boldsymbol{\theta}_j)\}}
#'        {\sum_{r=m}^{K_b}\exp\{u_{i_{(r)}}(\boldsymbol{\theta}_j)\}}.
#' }
#' MOLE and PICK probabilities are sums over compatible full rankings,
#' normalized over the observed block patterns by the C++ probability helper.
#'
#' \strong{Q-matrix for unfolding.} The Q-matrix serves a dual role:
#' entries of \eqn{\pm 1} activate both the slope \eqn{a_{id}} and the
#' location \eqn{\delta_{id}}; the sign of \eqn{q_{id}} determines the
#' sign of \eqn{\delta_{id}} (positive vs. negative side of dimension
#' \eqn{d}).
#'
#' @section Estimation Methods:
#'
#' \describe{
#'   \item{\strong{Stan} (\code{method = "stan"}):}{
#'     Full Bayesian inference via HMC. The GGUM likelihood involves
#'     summation terms that are handled in the Stan model block.}
#'   \item{\strong{iStEM} (\code{method = "iStEM"}):}{
#'     Improved Stochastic EM. Person sampling uses finite-grid block Gibbs;
#'     GGUM item parameters updated via constrained L-BFGS-B with
#'     monotonicity constraints on \eqn{\tau} and sign constraints on
#'     \eqn{\delta}.}
#' }
#'
#' @param data An \eqn{N \times B} forced-choice data matrix with ranking
#'   strings or 1-based pattern indices.
#' @param Q.matrix An optional \eqn{I \times D} matrix with entries
#'   \eqn{-1, 0, 1}. See \code{\link{fit.MGGUM}} for the sign convention.
#' @param block.items A list of length \eqn{B} giving global item indices
#'   in each block. Parsed from \code{data} if \code{NULL} and data contains
#'   ranking strings.
#' @param D Integer; number of latent dimensions. If \code{NULL}, inferred
#'   from \code{Q.matrix} when available, otherwise from the block structure.
#' @param fc.type Forced-choice response type: \code{"RANK"} (default),
#'   \code{"MOLE"}, or \code{"PICK"}.
#' @param method Estimation method: \code{"iStEM"} (default) or
#'   \code{"stan"}.
#' @param control.model A named list of model-level hyperparameters for the
#'   GGUM unfolding model in a forced-choice framework. Supported entries:
#'   \describe{
#'     \item{\code{a.mu}, \code{a.sigma}}{Prior location and scale for
#'           \eqn{\log a_{id}} (log-normal). Defaults: 0, 0.5.
#'           The log-normal prior ensures positivity of discrimination
#'           parameters in the unfolding model. The zero mean-log
#'           encourages slopes near 1 unless the data strongly indicate
#'           otherwise.}
#'     \item{\code{delta.pos.mu}, \code{delta.pos.sigma}}{Prior mean and
#'           SD for \eqn{\delta_{id}} when \eqn{q_{id} = 1} (positive-side
#'           item locations; normal). Defaults: 1, 0.5. In the FC context,
#'           these priors regularize the location of statements on the
#'           positive side of each latent dimension.}
#'     \item{\code{delta.neg.mu}, \code{delta.neg.sigma}}{Prior mean and
#'           SD for \eqn{\delta_{id}} when \eqn{q_{id} = -1} (negative-side
#'           item locations; normal). Defaults: -1, 0.5. These priors
#'           regularize statement locations on the negative side of each
#'           dimension.}
#'     \item{\code{tau.mu}, \code{tau.sigma}}{Location and scale for the
#'           log-normal prior on the positive threshold magnitude
#'           \eqn{-\tau_{i1}}. Defaults: 0, 0.5. In the binary FCGGUM,
#'           each statement has a single free threshold
#'           \eqn{\tau_{i1}<0}.}
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
#'           (NUTS; default: 0.95). The GGUM likelihood involves summation
#'           terms that can create challenging posterior geometries;
#'           consider increasing to 0.9--0.95 if divergences occur.}
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
#'     \item{\code{eps1}}{Geweke z-score convergence threshold for
#'           burn-in (default: 1.5).}
#'     \item{\code{eps2}}{Monte Carlo error tolerance for the final
#'           chain (default: 0.4).}
#'     \item{\code{frac1}, \code{frac2}}{Fractions for the Geweke
#'           diagnostic (defaults: 0.1, 0.5).}
#'     \item{\code{corr.optim.maxit}}{Maximum L-BFGS-B iterations for the
#'           constrained unit-diagonal correlation update (default: 50).}
#'     \item{\code{optim.maxit}}{Maximum L-BFGS-B iterations per item
#'           (default: 50). GGUM item parameters (\eqn{a}, \eqn{\delta},
#'           \eqn{\tau}) are optimized jointly under monotonicity and
#'           sign constraints. Increase if optimization warnings appear.}
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
#'           bounds on the positive gap used to construct
#'           \eqn{\tau_{i1}<0}. Defaults are 1e-4 and 6.}
#'   }
#'
#' @return An object of class \code{"FCGGUM"} containing:
#' \describe{
#'   \item{\code{npar}}{Number of free parameters.}
#'   \item{\code{method}}{\code{"stan"} or \code{"iStEM"}.}
#'   \item{\code{theta}}{List with \code{est}, \code{se}, \code{Rhat}
#'         (\eqn{N \times D}).}
#'   \item{\code{par}}{List with \code{est}, \code{se}, \code{Rhat},
#'         \code{free} (\eqn{I \times (2D + 2)}).}
#'   \item{\code{Corr}}{List with \code{est}, \code{se}, \code{Rhat}
#'         (\eqn{D \times D}).}
#'   \item{\code{block.items}, \code{patterns}, \code{patterns.total},
#'         \code{response}}{Block structure.}
#'   \item{\code{logLik}}{Marginal log-likelihood (class \code{"logLik"}).}
#' }
#'
#' @references
#' Roberts, J. S., Donoghue, J. R., & Laughlin, J. E. (2000). A general item
#'   response theory model for unfolding unidimensional polytomous responses.
#'   \emph{Applied Psychological Measurement}, 24(1), 3--32.
#'
#' @seealso
#' \code{\link{sim.data.FCGGUM}}, \code{\link{get.fit.index.FCGGUM}},
#' \code{\link{logLik.FCGGUM}}, \code{\link{fit.MGGUM}}
#'
#' @examples
#' sim <- sim.data.FCGGUM(N.person = 20, N.block = 3, I.block = 2,
#'                        D = 2, fc.type = "RANK")
#' fit <- fit.FCGGUM(sim$data, Q.matrix = sim$Q.matrix,
#'                   block.items = sim$block.items,
#'                   D = 2, fc.type = "RANK", method = "iStEM",
#'                   control.method = list(
#'                     vis = FALSE, seed = 123,
#'                     M = 2, B = 2, burnin.maxitr = 2,
#'                     maxitr = 3, eps1 = 10, eps2 = 10,
#'                     estimate.se = FALSE))
#' head(fit$theta$est)
#' gof <- get.fit.index(fit)
#' summary(gof)
#'
#' @export
fit.FCGGUM <- function(data, Q.matrix = NULL, block.items = NULL,
                       D = NULL, fc.type = NULL,
                       method = c("iStEM", "stan"),
                       control.model = NULL,
                       control.method = NULL) {

  call <- match.call()
  method <- match.arg(method)

  block.items <- resolve_block_items(block.items, data)
  fc.type     <- resolve_fc_type(fc.type, data, block.items)
  D           <- resolve_D(D, Q.matrix = Q.matrix, block.items = block.items)

  control.model  <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")

  if (method == "stan") {
    return(fit.FCGGUM.stan(
      data = data,
      Q.matrix = Q.matrix,
      block.items = block.items,
      D = D,
      fc.type = fc.type,
      control.model = control.model,
      control.method = control.method,
      .call = call
    ))
  }

  fit.FCGGUM.iStEM(
    data = data,
    Q.matrix = Q.matrix,
    block.items = block.items,
    D = D,
    fc.type = fc.type,
    control.model = control.model,
    control.method = control.method,
    .call = call
  )
}
