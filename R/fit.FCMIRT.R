#' Fit the Forced-Choice Multidimensional IRT (FCMIRT) Model
#'
#' @description
#' Fits a forced-choice variant of the multidimensional item response theory
#' model to comparative ranking data. In forced-choice (FC) formats,
#' respondents compare items within blocks rather than rating each item in
#' isolation. The FCMIRT model couples an item-level MIRT endorsement model
#' with a block-level ranking mechanism to recover latent trait estimates
#' from comparative responses.
#'
#' @section Model Specification:
#'
#' \strong{Item-level endorsement.}
#' Each statement \eqn{i = 1, \dots, I} is governed by a MIRT model
#' (1PL through 4PL). For a person with trait \eqn{\boldsymbol{\theta}_j},
#' the probability of endorsing statement \eqn{i} in isolation is:
#' \deqn{
#'   P_i(\boldsymbol{\theta}_j) =
#'   c_i + (d_i - c_i) \times
#'   \frac{1}{1 + \exp\bigl[-\bigl(\sum_{d=1}^{D} a_{id}\,\theta_{jd} - b_i\bigr)\bigr]}.
#' }
#'
#' \strong{Block-level ranking.}
#' Let block \eqn{b = 1, \dots, B} consist of items
#' \eqn{\{i_1, i_2, \dots, i_{K_b}\}} and define the implemented
#' statement utility as the logit of the item-level endorsement probability,
#' \deqn{
#'   u_i(\boldsymbol{\theta}_j) =
#'   \log\left[
#'     \frac{P_i(\boldsymbol{\theta}_j)}
#'          {1-P_i(\boldsymbol{\theta}_j)}
#'   \right].
#' }
#' The probability that an examinee
#' produces the ranking \eqn{i_{(1)} \succ i_{(2)} \succ \dots \succ i_{(K_b)}}
#' (read as "\eqn{i_{(1)}} is preferred over \eqn{i_{(2)}} over ...") is:
#' \deqn{
#'   P\bigl(i_{(1)} \succ \dots \succ i_{(K_b)} \mid \boldsymbol{\theta}_j\bigr) =
#'   \prod_{m=1}^{K_b-1}
#'   \frac{
#'     \exp\{u_{i_{(m)}}(\boldsymbol{\theta}_j)\}
#'   }{
#'     \sum_{r=m}^{K_b} \exp\{u_{i_{(r)}}(\boldsymbol{\theta}_j)\}
#'   }.
#' }
#' The Stan code evaluates the same choice kernel as
#' \eqn{\log P_i + \sum_{h \ne i}\log(1-P_h)} within each remaining set,
#' which is algebraically equivalent to a softmax over
#' \eqn{\mathrm{logit}(P_i)}.
#'
#' \strong{Partial rankings (MOLE / PICK).} For MOLE (most-least) data, only
#' the best and worst items in each block are identified; for PICK data,
#' only the best item. The probability of a partial ranking is obtained by
#' marginalizing (summing) the full-ranking probabilities over all
#' completions consistent with the partial constraint; the C++ probability
#' helper returns probabilities normalized over the observed RANK/MOLE/PICK
#' patterns for each block.
#'
#' \strong{Identification constraint (2PLM-RANK).} Block intercepts
#' \eqn{b_i} are constrained to sum to zero within each block:
#' \deqn{
#'   \sum_{i \in \text{block } b} b_i = 0, \qquad b = 1, \dots, B,
#' }
#' so that only \eqn{K_b - 1} difficulties are freely estimated per block.
#'
#' @section Estimation Methods:
#'
#' \describe{
#'   \item{\strong{Stan} (\code{method = "stan"}):}{
#'     Full Bayesian inference via HMC. The forced-choice likelihood is
#'     evaluated over the ranking pattern probabilities. By default,
#'     \code{a.sigma} is set to \code{1.0} (wider prior) to avoid spurious
#'     shrinkage of within-block slopes toward equality.}
#'   \item{\strong{iStEM} (\code{method = "iStEM"}):}{
#'     Improved Stochastic EM with finite-grid block Gibbs person-sampling
#'     and item optimization (L-BFGS-B). Block-level identification
#'     constraints are enforced during each item-parameter update.}
#' }
#'
#' @param data An \eqn{N \times B} forced-choice data object. Each entry is
#'   either a ranking string (e.g., \code{"2>1>3"}) or a 1-based pattern
#'   index into the block's ranking patterns. If pattern indices are used,
#'   \code{block.items} must be supplied.
#' @param model MIRT model type: \code{"m1pl"}, \code{"m2pl"} (default),
#'   \code{"m3pl"}, or \code{"m4pl"}.
#' @param Q.matrix An optional \eqn{I \times D} binary (0/1) matrix
#'   indicating active item loadings.
#' @param block.items A list of length \eqn{B} where each element is an
#'   integer vector of global item indices belonging to that block. If
#'   \code{NULL} and \code{data} contains ranking strings, automatically
#'   parsed from the data.
#' @param D Integer; number of latent dimensions. If \code{NULL}, inferred
#'   from \code{Q.matrix}; if both are \code{NULL}, an error is raised.
#' @param fc.type Character vector (length 1 or \eqn{B}): \code{"RANK"}
#'   (full ranking, default), \code{"MOLE"} (most-least), or \code{"PICK"}
#'   (best only). Recycled to all blocks if scalar.
#' @param method Estimation method: \code{"iStEM"} (default) or
#'   \code{"stan"}.
#' @param control.model A named list of model-level hyperparameters.
#'   Supported entries:
#'   \describe{
#'     \item{\code{a.mu}, \code{a.sigma}}{Prior location and scale for
#'           \eqn{\log a_{id}} (log-normal). Defaults: 0, 1.0.
#'           \strong{Important:} The Stan default for \code{a.sigma}
#'           is deliberately set to 1.0 (wider than the MIRT default of
#'           0.25) because forced-choice blocks identify relative
#'           utilities; tight common slope priors can spuriously
#'           shrink within-block slopes toward equality, degrading
#'           trait recovery. For iStEM, the prior defaults match the
#'           standard MIRT values (\code{a.mu = 0.25}, \code{a.sigma = 0.25})
#'           unless overridden.}
#'     \item{\code{b.mu}, \code{b.sigma}}{Prior mean and SD for \eqn{b_i}
#'           (normal). Defaults: 0, 1. Under the block zero-sum constraint
#'           \eqn{\sum_{i \in \text{block } b} b_i = 0}, these priors
#'           regularize the free difficulty parameters.}
#'     \item{\code{c.mu}, \code{c.sigma}}{Uniform support bounds
#'           \eqn{[c_{\min}, c_{\max}]} for lower-asymptote parameters
#'           \eqn{c_i}. Defaults: 0, 0.35. Only used for 3PL/4PL models.}
#'     \item{\code{d.mu}, \code{d.sigma}}{Uniform support bounds
#'           \eqn{[d_{\min}, d_{\max}]} for upper-asymptote parameters
#'           \eqn{d_i}. Defaults: 0.65, 1. Only used for 4PL models.}
#'     \item{\code{theta.mu}}{Prior mean vector for
#'           \eqn{\boldsymbol{\theta}_j}. Default: \code{rep(0, D)}.}
#'     \item{\code{L}}{Theta grid size per dimension for marginal
#'           log-likelihood computation and iStEM block Gibbs sampling.
#'           Default adapts to \eqn{D}
#'           (e.g., 61 for D = 1, 31 for D = 2, 15 for D = 3).}
#'     \item{\code{theta.lower}, \code{theta.upper}}{Bounds for the
#'           theta grid used by marginal log-likelihood computation and
#'           iStEM block Gibbs sampling. Defaults: -6, 6.}
#'     \item{\code{use.prior}}{Logical (iStEM only). If \code{FALSE}, the
#'           item-prior penalty term is omitted from the item update,
#'           recovering an approximate maximum-likelihood item step.
#'           Default is \code{TRUE}.}
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
#'           (NUTS; default: 0.95). Forced-choice models often have
#'           challenging posterior geometries; consider increasing to
#'           0.9 or higher if divergent transitions appear.}
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
#'           (default: 50).}
#'     \item{\code{fix.corr}}{Logical; fix correlations to identity
#'           (default: \code{FALSE}).}
#'     \item{\code{estimate.se}}{Logical; compute standard errors from
#'           final MC chain (default: \code{TRUE}).}
#'     \item{\code{a.lower}, \code{a.upper}}{Bounds on \eqn{a_{id}}
#'           during optimization (defaults: 1e-4, 6).}
#'     \item{\code{b.lower}, \code{b.upper}}{Bounds on \eqn{b_i}
#'           during optimization (defaults: -6, 6).}
#'     \item{\code{c.lower}, \code{c.upper}}{Bounds on \eqn{c_i} for
#'           3PL/4PL models. Defaults from prior support.}
#'     \item{\code{d.lower}, \code{d.upper}}{Bounds on \eqn{d_i} for
#'           4PL models. Defaults from prior support.}
#'   }
#'
#' @return An object of class \code{"FCMIRT"} containing:
#' \describe{
#'   \item{\code{npar}}{Number of free parameters.}
#'   \item{\code{method}}{\code{"stan"} or \code{"iStEM"}.}
#'   \item{\code{theta}}{List with \code{est}, \code{se}, \code{Rhat}
#'         (\eqn{N \times D}); person trait estimates.}
#'   \item{\code{par}}{List with \code{est}, \code{se}, \code{Rhat},
#'         \code{free} (\eqn{I \times (D+3)}); item parameters
#'         \code{a1..aD, b, c, d}.}
#'   \item{\code{Corr}}{List with \code{est}, \code{se}, \code{Rhat}
#'         (\eqn{D \times D}).}
#'   \item{\code{block.items}, \code{patterns}, \code{patterns.total},
#'         \code{response}}{Block structure and response data.}
#'   \item{\code{logLik}}{Marginal log-likelihood (class \code{"logLik"}).}
#' }
#'
#' @references
#' Brown, A., & Maydeu-Olivares, A. (2011). Item response modeling of
#'   forced-choice questionnaires. \emph{Educational and Psychological
#'   Measurement}, 71(3), 460--502. \doi{10.1177/0013164410375112}
#'
#' Luce, R. D. (1959). \emph{Individual choice behavior: A theoretical
#'   analysis}. Wiley.
#'
#' Plackett, R. L. (1975). The analysis of permutations. \emph{Journal of
#'   the Royal Statistical Society: Series C}, 24(2), 193--202.
#'
#' @seealso
#' \code{\link{sim.data.FCMIRT}}, \code{\link{get.fit.index.FCMIRT}},
#' \code{\link{logLik.FCMIRT}}, \code{\link{fit.MIRT}}
#'
#' @examples
#' # Simulate forced-choice data
#' sim <- sim.data.FCMIRT(N.person = 20, N.block = 3, I.block = 2,
#'                        D = 2, model = "m2pl", fc.type = "RANK")
#'
#' # Fit via iStEM
#' fit <- fit.FCMIRT(sim$data, model = "m2pl",
#'                   Q.matrix = sim$Q.matrix,
#'                   block.items = sim$block.items,
#'                   D = 2, fc.type = "RANK",
#'                   method = "iStEM",
#'                   control.method = list(
#'                     vis = FALSE, seed = 123,
#'                     M = 2, B = 2, burnin.maxitr = 2,
#'                     maxitr = 3, eps1 = 10, eps2 = 10,
#'                     estimate.se = FALSE))
#'
#' # Examine trait recovery
#' cor(fit$theta$est, sim$theta)
#'
#' # Fit indices (nominal binary expansion)
#' gof <- get.fit.index(fit)
#' summary(gof)
#'
#' @export
fit.FCMIRT <- function(data, model = "2PL", Q.matrix = NULL,
                       block.items = NULL, D = NULL,
                       fc.type = NULL,
                       method = c("iStEM", "stan"),
                       control.model = NULL,
                       control.method = NULL) {

  call <- match.call()
  method <- match.arg(method)
  model <- resolve_model_type(model)

  # Auto-detect block.items and fc.type from character ranking data
  block.items <- resolve_block_items(block.items, data)
  fc.type     <- resolve_fc_type(fc.type, data, block.items)
  D           <- resolve_D(D, Q.matrix = Q.matrix, block.items = block.items)

  control.model  <- fc_as_control_list(control.model, "control.model")
  control.method <- fc_as_control_list(control.method, "control.method")

  if (method == "stan") {
    return(fit.FCMIRT.stan(
      data = data,
      model = model,
      Q.matrix = Q.matrix,
      block.items = block.items,
      D = D,
      fc.type = fc.type,
      control.model = control.model,
      control.method = control.method,
      .call = call
    ))
  }

  fit.FCMIRT.iStEM(
    data = data,
    model = model,
    Q.matrix = Q.matrix,
    block.items = block.items,
    D = D,
    fc.type = fc.type,
    control.model = control.model,
    control.method = control.method,
    .call = call
  )
}
