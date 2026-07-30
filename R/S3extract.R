#' @title S3 Methods: extract
#'
#' @description
#' A generic S3 extractor function designed to retrieve internal components
#' from fitted model objects produced by the \pkg{ForceChoice} package. This
#' function provides a consistent interface across all eight model classes,
#' allowing users to access estimated parameters, fit statistics, latent
#' trait estimates, convergence diagnostics, and model configuration data
#' without directly manipulating the internal list structure.
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
#' @param what A character string specifying the name of the component to
#'   extract. Valid choices depend on the class of \code{object}. See Details
#'   section for full listings.
#' @param ... Additional arguments passed to methods (currently ignored).
#'
#' @return The requested component. Return type varies depending on
#'   \code{what} and the class of \code{object}. If an invalid \code{what}
#'   is provided, an informative error is thrown listing valid options.
#'
#' @details
#' This function supports extraction from the eight ForceChoice model classes.
#' Below are the available components for each:
#'
#' \describe{
#'   \item{\strong{\code{MIRT}, \code{MGPCM}, \code{MGGUM}}}{
#'     Traditional (single-stimulus) models. Available components:
#'     \describe{
#'       \item{\code{par}}{List (\code{est}, \code{se}, \code{Rhat}, \code{free})
#'         of item parameter arrays.}
#'       \item{\code{theta}}{List (\code{est}, \code{se}, \code{Rhat}) of
#'         \eqn{N \times D} person parameter matrices.}
#'       \item{\code{Corr}}{List (\code{est}, \code{se}, \code{Rhat}) of
#'         \eqn{D \times D} inter-trait correlation matrices.}
#'       \item{\code{logLik}}{Marginal log-likelihood (class \code{"logLik"}).}
#'       \item{\code{npar}}{Number of free parameters.}
#'       \item{\code{method}}{Estimation method (\code{"stan"} or \code{"iStEM"}).}
#'       \item{\code{Q.matrix}}{The \eqn{I \times D} Q-matrix.}
#'       \item{\code{length.poly}}{Vector of category counts per item
#'         (MGPCM, MGGUM only).}
#'       \item{\code{stan.obj}}{The \code{stanfit} object (Stan only).}
#'       \item{\code{MCMC.obj}}{The MCMC extract list (Stan only).}
#'       \item{\code{iStEM}}{List of iStEM diagnostics (iStEM only).}
#'       \item{\code{call}}{The matched call.}
#'       \item{\code{arguments}}{List of arguments used for fitting.}
#'     }
#'   }
#'
#'   \item{\strong{\code{FCMIRT}, \code{FCGGUM}}}{
#'     Forced-choice models. In addition to the traditional-model components
#'     listed above, the following are available:
#'     \describe{
#'       \item{\code{block.items}}{List of item indices in each block.}
#'       \item{\code{response}}{The \eqn{N \times B} forced-choice response matrix.}
#'       \item{\code{patterns}}{List of permissible ranking patterns per block.}
#'       \item{\code{patterns.total}}{List of all possible ranking patterns per block.}
#'       \item{\code{fc.type}}{Character vector of forced-choice formats
#'         (\code{"RANK"}, \code{"MOLE"}, \code{"PICK"}).}
#'     }
#'   }
#'
#'   \item{\strong{\code{FCDCM}}}{
#'     Forced-Choice DCM. Available components:
#'     \describe{
#'       \item{\code{par}}{List (\code{est}, \code{se}, \code{Rhat}, \code{free})
#'         of block \eqn{\eta} parameter matrices (\eqn{B \times 2}).}
#'       \item{\code{delta}}{List (\code{est}, \code{se}, \code{Rhat}) of
#'         higher-order \eqn{\delta} parameters (\eqn{D \times 2}).}
#'       \item{\code{theta}}{List (\code{est}, \code{se}, \code{Rhat}) of
#'         \eqn{N \times 1} higher-order trait estimates.}
#'       \item{\code{alpha}}{Posterior mean attribute profile
#'         (\eqn{N \times D}).}
#'       \item{\code{alpha.patterns}}{Full enumeration of \eqn{2^D}
#'         attribute mastery patterns.}
#'       \item{\code{zeta.patterns}}{Condensation outputs for each pattern.}
#'       \item{\code{dcm.type}}{DCM condensation rule (\code{"DINA"} or
#'         \code{"DINO"}).}
#'       \item{\code{response}}{The \eqn{N \times B} binary block response matrix.}
#'       \item{\code{block.items}}{List of two-statement blocks.}
#'       \item{\code{patterns}}{List of binary response patterns per block.}
#'       \item{\code{logLik}}{Marginal log-likelihood (class \code{"logLik"}).}
#'       \item{\code{npar}}{Number of free parameters.}
#'       \item{\code{method}}{Estimation method.}
#'       \item{\code{Q.matrix}}{The \eqn{I \times D} Q-matrix.}
#'       \item{\code{stan.obj}, \code{MCMC.obj}, \code{iStEM}, \code{EM}}{
#'         Estimation backend objects.}
#'       \item{\code{call}}{The matched call.}
#'       \item{\code{arguments}}{List of arguments used for fitting.}
#'     }
#'   }
#'
#'   \item{\strong{\code{FCGDINA}}}{
#'     Forced-Choice GDINA. Available components:
#'     \describe{
#'       \item{\code{delta}}{List (\code{est}, \code{se}, \code{Rhat}) of
#'         item-level CDM delta parameters.}
#'       \item{\code{alpha}}{Posterior attribute summaries, including
#'         posterior class probabilities when available.}
#'       \item{\code{alpha.patterns}}{Full enumeration of \eqn{2^D}
#'         attribute mastery patterns.}
#'       \item{\code{design.matrix.list}}{Per-item CDM design matrices.}
#'       \item{\code{response}, \code{block.items}, \code{patterns},
#'         \code{patterns.total}, \code{fc.type}}{Forced-choice data
#'         structures.}
#'     }
#'   }
#'
#'   \item{\strong{\code{TIRT}}}{
#'     Thurstonian IRT. Available components:
#'     \describe{
#'       \item{\code{par}}{List (\code{est}, \code{se}, \code{Rhat}, \code{free})
#'         of statement-level parameter matrices (\eqn{I \times (D+1)}).}
#'       \item{\code{theta}}{List (\code{est}, \code{se}, \code{Rhat}) of
#'         \eqn{N \times D} person parameter matrices.}
#'       \item{\code{gamma}}{List (\code{est}, \code{se}, \code{Rhat},
#'         \code{free}) of pairwise \eqn{\gamma} parameters.}
#'       \item{\code{gamma.matrix}}{List (\code{est}, \code{se}, \code{Rhat})
#'         of \eqn{I \times I} skew-symmetric \eqn{\gamma} matrices.}
#'       \item{\code{Corr}}{List (\code{est}, \code{se}, \code{Rhat}) of
#'         \eqn{D \times D} inter-trait correlation matrices.}
#'       \item{\code{pairs.matrix}}{Matrix of pairwise comparison indices.}
#'       \item{\code{pairs.value}}{Person-specific pair data (MOLE/PICK).}
#'       \item{\code{response}}{The \eqn{N \times I_{pairs}} pairwise
#'         response matrix.}
#'       \item{\code{block.items}}{List of item indices in each block.}
#'       \item{\code{fc.type}}{Character vector of forced-choice formats.}
#'       \item{\code{Q.matrix}}{The statement-level Q-matrix.}
#'       \item{\code{logLik}}{Marginal log-likelihood (class \code{"logLik"}).}
#'       \item{\code{npar}}{Number of free parameters.}
#'       \item{\code{method}}{Estimation method.}
#'       \item{\code{stan.obj}, \code{MCMC.obj}, \code{iStEM}}{Estimation
#'         backend objects.}
#'       \item{\code{call}}{The matched call.}
#'       \item{\code{arguments}}{List of arguments used for fitting.}
#'     }
#'   }
#' }
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
#'
#' extract(fit, "par")       # item parameter estimates
#' extract(fit, "theta")     # person trait estimates
#' extract(fit, "Corr")      # factor correlation matrix
#' extract(fit, "npar")      # number of free parameters
#' extract(fit, "logLik")    # marginal log-likelihood
#' extract(fit, "iStEM")     # iStEM convergence diagnostics
#'
#' @name extract
NULL

#' @export
extract <- function(object, what, ...) {
  UseMethod("extract")
}


# ===========================================================================
# extract.MIRT
# ===========================================================================

#' @describeIn extract Extract components from \code{MIRT} objects.
#' @method extract MIRT
#' @export
extract.MIRT <- function(object, what, ...) {
  choices <- c("par", "theta", "Corr", "logLik", "npar", "method",
               "Q.matrix", "stan.obj", "MCMC.obj", "iStEM",
               "call", "arguments")
  what <- match.arg(what, choices)
  switch(what,
    par       = object$par,
    theta     = object$theta,
    Corr      = object$Corr,
    logLik    = object$logLik,
    npar      = object$npar,
    method    = object$method,
    Q.matrix  = object$Q.matrix,
    stan.obj  = object$stan.obj,
    MCMC.obj  = object$MCMC.obj,
    iStEM     = object$iStEM,
    call      = object$call,
    arguments = object$arguments
  )
}

# ===========================================================================
# extract.MGPCM
# ===========================================================================

#' @describeIn extract Extract components from \code{MGPCM} objects.
#' @method extract MGPCM
#' @export
extract.MGPCM <- function(object, what, ...) {
  choices <- c("par", "theta", "Corr", "logLik", "npar", "method",
               "Q.matrix", "length.poly",
               "stan.obj", "MCMC.obj", "iStEM",
               "call", "arguments")
  what <- match.arg(what, choices)
  switch(what,
    par         = object$par,
    theta       = object$theta,
    Corr        = object$Corr,
    logLik      = object$logLik,
    npar        = object$npar,
    method      = object$method,
    Q.matrix    = object$Q.matrix,
    length.poly = object$length.poly,
    stan.obj    = object$stan.obj,
    MCMC.obj    = object$MCMC.obj,
    iStEM       = object$iStEM,
    call        = object$call,
    arguments   = object$arguments
  )
}

# ===========================================================================
# extract.MGGUM
# ===========================================================================

#' @describeIn extract Extract components from \code{MGGUM} objects.
#' @method extract MGGUM
#' @export
extract.MGGUM <- function(object, what, ...) {
  choices <- c("par", "theta", "Corr", "logLik", "npar", "method",
               "Q.matrix", "length.poly",
               "stan.obj", "MCMC.obj", "iStEM",
               "call", "arguments")
  what <- match.arg(what, choices)
  switch(what,
    par         = object$par,
    theta       = object$theta,
    Corr        = object$Corr,
    logLik      = object$logLik,
    npar        = object$npar,
    method      = object$method,
    Q.matrix    = object$Q.matrix,
    length.poly = object$length.poly,
    stan.obj    = object$stan.obj,
    MCMC.obj    = object$MCMC.obj,
    iStEM       = object$iStEM,
    call        = object$call,
    arguments   = object$arguments
  )
}

# ===========================================================================
# extract.FCMIRT
# ===========================================================================

#' @describeIn extract Extract components from \code{FCMIRT} objects.
#' @method extract FCMIRT
#' @export
extract.FCMIRT <- function(object, what, ...) {
  choices <- c("par", "theta", "Corr", "logLik", "npar", "method",
               "Q.matrix", "block.items", "response",
               "patterns", "patterns.total", "fc.type",
               "stan.obj", "MCMC.obj", "iStEM",
               "call", "arguments")
  what <- match.arg(what, choices)
  switch(what,
    par            = object$par,
    theta          = object$theta,
    Corr           = object$Corr,
    logLik         = object$logLik,
    npar           = object$npar,
    method         = object$method,
    Q.matrix       = object$Q.matrix,
    block.items    = object$block.items,
    response       = object$response,
    patterns       = object$patterns,
    patterns.total = object$patterns.total,
    fc.type        = object$fc.type,
    stan.obj       = object$stan.obj,
    MCMC.obj       = object$MCMC.obj,
    iStEM          = object$iStEM,
    call           = object$call,
    arguments      = object$arguments
  )
}

# ===========================================================================
# extract.FCDCM
# ===========================================================================

#' @describeIn extract Extract components from \code{FCDCM} objects.
#' @method extract FCDCM
#' @export
extract.FCDCM <- function(object, what, ...) {
  choices <- c("par", "delta", "theta", "alpha", "class.post",
               "alpha.patterns", "zeta.patterns", "dcm.type",
               "logLik", "npar", "method",
               "Q.matrix", "block.items", "response", "patterns",
               "stan.obj", "MCMC.obj", "iStEM",
               "call", "arguments")
  what <- match.arg(what, choices)
  switch(what,
    par             = object$par,
    delta           = object$delta,
    theta           = object$theta,
    alpha           = object$alpha,
    class.post      = object$class.post,
    alpha.patterns  = object$alpha.patterns,
    zeta.patterns   = object$zeta.patterns,
    dcm.type        = object$dcm.type,
    logLik          = object$logLik,
    npar            = object$npar,
    method          = object$method,
    Q.matrix        = object$Q.matrix,
    block.items     = object$block.items,
    response        = object$response,
    patterns        = object$patterns,
    stan.obj        = object$stan.obj,
    MCMC.obj        = object$MCMC.obj,
    iStEM           = object$iStEM,
    call            = object$call,
    arguments       = object$arguments
  )
}

# ===========================================================================
# extract.FCGDINA
# ===========================================================================

#' @describeIn extract Extract components from \code{FCGDINA} objects.
#' @method extract FCGDINA
#' @export
extract.FCGDINA <- function(object, what, ...) {
  choices <- c("delta", "alpha", "alpha.patterns", "design.matrix.list",
               "logLik", "npar", "method", "Q.matrix", "block.items",
               "response", "patterns", "patterns.total", "fc.type",
               "stan.obj", "MCMC.obj", "iStEM", "EM", "call", "arguments")
  what <- match.arg(what, choices)
  switch(what,
    delta             = object$delta,
    alpha             = object$alpha,
    alpha.patterns    = object$alpha.patterns,
    design.matrix.list = object$design.matrix.list,
    logLik            = object$logLik,
    npar              = object$npar,
    method            = object$method,
    Q.matrix          = object$Q.matrix,
    block.items       = object$block.items,
    response          = object$response,
    patterns          = object$patterns,
    patterns.total    = object$patterns.total,
    fc.type           = object$fc.type,
    stan.obj          = object$stan.obj,
    MCMC.obj          = object$MCMC.obj,
    iStEM             = object$iStEM,
    EM                = object$EM,
    call              = object$call,
    arguments         = object$arguments
  )
}

# ===========================================================================
# extract.FCGGUM
# ===========================================================================

#' @describeIn extract Extract components from \code{FCGGUM} objects.
#' @method extract FCGGUM
#' @export
extract.FCGGUM <- function(object, what, ...) {
  choices <- c("par", "theta", "Corr", "logLik", "npar", "method",
               "Q.matrix", "block.items", "response", "length.poly",
               "patterns", "patterns.total", "fc.type",
               "stan.obj", "MCMC.obj", "iStEM",
               "call", "arguments")
  what <- match.arg(what, choices)
  switch(what,
    par            = object$par,
    theta          = object$theta,
    Corr           = object$Corr,
    logLik         = object$logLik,
    npar           = object$npar,
    method         = object$method,
    Q.matrix       = object$Q.matrix,
    block.items    = object$block.items,
    response       = object$response,
    length.poly    = object$length.poly,
    patterns       = object$patterns,
    patterns.total = object$patterns.total,
    fc.type        = object$fc.type,
    stan.obj       = object$stan.obj,
    MCMC.obj       = object$MCMC.obj,
    iStEM          = object$iStEM,
    call           = object$call,
    arguments      = object$arguments
  )
}

# ===========================================================================
# extract.TIRT
# ===========================================================================

#' @describeIn extract Extract components from \code{TIRT} objects.
#' @method extract TIRT
#' @export
extract.TIRT <- function(object, what, ...) {
  choices <- c("par", "theta", "gamma", "gamma.matrix", "Corr",
               "logLik", "npar", "method",
               "Q.matrix", "pairs.matrix", "pairs.value",
               "block.items", "response", "fc.type",
               "stan.obj", "MCMC.obj", "iStEM",
               "call", "arguments")
  what <- match.arg(what, choices)
  switch(what,
    par          = object$par,
    theta        = object$theta,
    gamma        = object$gamma,
    gamma.matrix = object$gamma.matrix,
    Corr         = object$Corr,
    logLik       = object$logLik,
    npar         = object$npar,
    method       = object$method,
    Q.matrix     = object$Q.matrix,
    pairs.matrix = object$pairs.matrix,
    pairs.value  = object$pairs.value,
    block.items  = object$block.items,
    response     = object$response,
    fc.type      = object$fc.type,
    stan.obj     = object$stan.obj,
    MCMC.obj     = object$MCMC.obj,
    iStEM        = object$iStEM,
    call         = object$call,
    arguments    = object$arguments
  )
}
