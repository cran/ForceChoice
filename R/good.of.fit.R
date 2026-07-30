#' Goodness-of-Fit Indices for Latent Probability Models
#'
#' @description
#' Computes a comprehensive suite of model-fit diagnostics for probability
#' models used in IRT, CDM, latent class, and forced-choice applications.
#' Indices include likelihood-based information criteria, limited-information
#' goodness-of-fit statistics (M2, RMSEA, SRMSR), incremental/comparative fit
#' indices (CFI, TLI), pseudo-\eqn{R^2} measures, residual diagnostics,
#' local-dependence statistics (Yen's Q3), and posterior classification
#' diagnostics.
#'
#' @section Likelihood and Information Criteria:
#'
#' Let \eqn{\ell = \log L(\hat{\boldsymbol{\psi}})} be the maximized
#' log-likelihood, \eqn{N} the sample size, and \eqn{p} the number of free
#' parameters (\code{npar}).
#'
#' \strong{Deviance:}
#' \deqn{D = -2\ell}
#'
#' \strong{AIC} (Akaike, 1974):
#' \deqn{\text{AIC} = D + 2p}
#'
#' \strong{AICc} (Hurvich & Tsai, 1989), small-sample correction:
#' \deqn{\text{AICc} = \text{AIC} + \frac{2p(p+1)}{N - p - 1},
#'       \quad N > p + 1}
#'
#' \strong{BIC} (Schwarz, 1978):
#' \deqn{\text{BIC} = D + p \log N}
#'
#' \strong{CAIC} (Bozdogan, 1987), consistent AIC:
#' \deqn{\text{CAIC} = D + p(\log N + 1)}
#'
#' \strong{SABIC} (Sclove, 1987), sample-size adjusted BIC:
#' \deqn{\text{SABIC} = D + p \log\bigl(\frac{N + 2}{24}\bigr)}
#'
#' \strong{HQIC} (Hannan & Quinn, 1979):
#' \deqn{\text{HQIC} = D + 2p \log(\log N), \quad N > e}
#'
#' @section Limited-Information M2 Family:
#'
#' The limited-information framework (Maydeu-Olivares & Joe, 2005, 2006) uses
#' first- and second-order marginal moments rather than the full
#' \eqn{I}-way contingency table. This is essential for sparse high-dimensional
#' categorical responses where full-information tests break down.
#'
#' Let \eqn{\mathbf{m}(\boldsymbol{\psi})} be the vector of model-implied
#' univariate and bivariate moments (dimension \eqn{M}), and let
#' \eqn{\hat{\boldsymbol{\Xi}}_2} be their asymptotic covariance matrix
#' under the model.
#'
#' \strong{M2 statistic:}
#' \deqn{
#'   M_2 = N\,
#'   \mathbf{e}'\hat{\mathbf{C}}_2\mathbf{e},
#' }
#' where \eqn{\mathbf{e} = \mathbf{p} - \boldsymbol{\pi}(\hat{\boldsymbol{\psi}})}
#' is the vector of marginal residuals, and
#' \eqn{\hat{\mathbf{C}}_2 = \boldsymbol{\Delta}_c(\boldsymbol{\Delta}_c'
#' \hat{\boldsymbol{\Xi}}_2 \boldsymbol{\Delta}_c)^{-1}\boldsymbol{\Delta}_c'}
#' with \eqn{\boldsymbol{\Delta}_c} the orthogonal complement of the Jacobian
#' matrix \eqn{\boldsymbol{\Delta} = \partial\boldsymbol{\pi}/\partial\boldsymbol{\psi}'}.
#'
#' Under correct model specification, \eqn{M_2 \sim \chi^2_{M - p}}
#' asymptotically.
#'
#' \strong{RMSEA (Root Mean Square Error of Approximation):}
#' \deqn{
#'   \text{RMSEA} = \sqrt{\max\!\left(0,\; \frac{M_2 - df}{N \cdot df}\right)},
#'   \qquad df = M - p
#' }
#' Confidence intervals are obtained via non-central \eqn{\chi^2} inversion.
#' The close-fit test evaluates \eqn{H_0: \text{RMSEA} \le 0.05}.
#'
#' \strong{SRMSR (Standardized Root Mean Square Residual):}
#' \deqn{
#'   \text{SRMSR} = \sqrt{
#'     \frac{1}{J} \sum_{i < j} (r_{ij}^{\text{obs}} - r_{ij}^{\text{model}})^2
#'   },
#' }
#' where \eqn{r_{ij}} are the observed and model-implied correlation residuals,
#' summed over \eqn{J} unique off-diagonal pairs.
#'
#' \strong{McDonald's NCI (Non-Centrality Index):}
#' \deqn{
#'   \text{NCI} = \exp\!\left(-\frac{\max(M_2 - df, 0)}{2N}\right)
#' }
#'
#' @section Incremental / Comparative Fit Indices:
#'
#' \strong{CFI (Comparative Fit Index):}
#' \deqn{
#'   \text{CFI} = 1 - \frac{\max(M_2 - df, 0)}
#'                           {\max(M_2^{\text{null}} - df^{\text{null}}, M_2 - df, 0)},
#' }
#' where the null (independence) model has \eqn{M_2^{\text{null}}} with
#' \eqn{df^{\text{null}}} degrees of freedom.
#'
#' \strong{TLI (Tucker--Lewis Index, a.k.a. NNFI):}
#' \deqn{
#'   \text{TLI} = \frac{M_2^{\text{null}} / df^{\text{null}} - M_2 / df}
#'                     {M_2^{\text{null}} / df^{\text{null}} - 1}
#' }
#'
#' \strong{IFI (Incremental Fit Index, a.k.a. BFI):}
#' \deqn{
#'   \text{IFI} = \frac{M_2^{\text{null}} - M_2}{M_2^{\text{null}} - df}
#' }
#'
#' @section Pseudo-\eqn{R^2} Measures:
#'
#' \strong{McFadden's \eqn{R^2}:}
#' \deqn{R^2_{\text{McF}} = 1 - \frac{\ell}{\ell_0}}
#'
#' \strong{McFadden's adjusted \eqn{R^2}:}
#' \deqn{R^2_{\text{McF,adj}} = 1 - \frac{\ell - p}{\ell_0}}
#'
#' \strong{Cox--Snell \eqn{R^2}:}
#' \deqn{R^2_{\text{CS}} = 1 - \exp\!\left(\frac{2}{N}(\ell_0 - \ell)\right)}
#'
#' \strong{Nagelkerke \eqn{R^2}:}
#' \deqn{R^2_{\text{N}} = \frac{R^2_{\text{CS}}}
#'                                 {1 - \exp(2\ell_0 / N)}}
#'
#' \strong{Aldrich--Nelson \eqn{R^2}:}
#' \deqn{R^2_{\text{AN}} = \frac{G^2}{G^2 + N},
#'       \quad G^2 = D_0 - D}
#'
#' \strong{Veall--Zimmermann \eqn{R^2}:}
#' \deqn{R^2_{\text{VZ}} = R^2_{\text{AN}} \times
#'       \frac{-2\ell_0 + N}{-2\ell_0}}
#'
#' where \eqn{\ell_0} is the null-model (independence) log-likelihood and
#' \eqn{D_0} the corresponding deviance.
#'
#' @section Yen's Q3 Local Dependence:
#'
#' For each pair of items \eqn{(i, k)}, the residual correlation is
#' \deqn{
#'   Q_{3,ik} = \text{Corr}(r_{ji}, r_{jk}),
#'   \qquad r_{ji} = Y_{ji} - E[Y_{ji} \mid \hat{\boldsymbol{\theta}}_j],
#' }
#' where residuals are computed from posterior expected scores. Adjusted
#' Q3 values subtract the mean Q3 across all pairs to center the
#' distribution (Christensen, Makransky, & Horton, 2017).
#'
#' @section Heuristic Guidelines:
#'
#' The \code{print.summary.good.of.fit()} method reports descriptive
#' recommended ranges based on common conventions. These are heuristics,
#' not decision rules; fit cutoffs depend on model family, dimensionality,
#' category sparsity, local dependence, and sample size.
#' \itemize{
#'   \item \strong{RMSEA}: \eqn{\le 0.05} close fit; \eqn{\le 0.08} reasonable;
#'         \eqn{\le 0.10} mediocre (Browne & Cudeck, 1993)
#'   \item \strong{SRMSR}: \eqn{\le 0.08} acceptable (Hu & Bentler, 1999)
#'   \item \strong{CFI / TLI}: \eqn{\ge 0.95} good; \eqn{\ge 0.90} acceptable
#' }
#'
#' @param par.vec Numeric vector of free model parameters.
#' @param loglik.fun Function returning the fitted model log-likelihood.
#' @param prob.fun Function returning conditional response probabilities at
#'   the latent support points.
#' @param response An \eqn{N \times I} response matrix (binary, polytomous,
#'   or nominal coding depending on \code{response.type}).
#' @param npar Integer; number of free parameters in the model.
#' @param pi Numeric vector of latent support weights (summing to 1).
#' @param pi.fun Optional function returning latent support weights given
#'   \code{par.vec} and \code{...}. If \code{NULL}, \code{pi} is used as-is.
#' @param alpha Tail probability for the RMSEA confidence interval.
#'   Default \code{0.05} yields a 90% CI (two-sided \eqn{2\alpha}).
#'   Use \code{0.025} for 95% or \code{0.005} for 99%.
#' @param response.type Character string specifying the response coding:
#'   \code{"auto"} (default; detects from probability dimensions),
#'   \code{"binary"} (each column is a binary indicator),
#'   \code{"polytomous"} (stacked category probabilities, requires
#'   \code{response.K}), or \code{"nominal"} (expanded binary indicators
#'   with mutually exclusive categories within groups, requires
#'   \code{nominal.groups}).
#' @param response.K Integer vector of length \eqn{I} giving the number of
#'   categories for each response column. Required when
#'   \code{response.type = "polytomous"}.
#' @param bivariate.groups Optional integer vector of length \eqn{I};
#'   response column pairs within the same group are excluded from bivariate
#'   limited-information moments. Useful for forced-choice or testlet
#'   designs where within-block comparisons are structural zeros.
#' @param nominal.groups Integer vector of length matching the number of
#'   expanded binary indicators for \code{response.type = "nominal"}.
#'   Indicators sharing the same group value belong to the same original
#'   nominal item/block.
#' @param exclusive.groups Deprecated alias for \code{nominal.groups}.
#' @param n.boot Integer; number of bootstrap replicates for bootstrapped
#'   diagnostics. Default \code{1000}; set to \code{0} to skip.
#' @param ... Additional arguments passed to \code{loglik.fun},
#'   \code{prob.fun}, and \code{pi.fun}.
#'
#' @return An object of class \code{"good.of.fit"} (a list). Key components:
#' \describe{
#'   \item{\code{LogLik}, \code{deviance}}{Log-likelihood and deviance.}
#'   \item{\code{npar}, \code{N}, \code{nitems}}{Parameter count, sample
#'         size, number of items/indicators.}
#'   \item{\code{AIC}, \code{AICc}, \code{BIC}, \code{CAIC}, \code{SABIC},
#'         \code{HQIC}}{Information criteria.}
#'   \item{\code{M2.value}, \code{df}, \code{p.value}}{M2 statistic, degrees
#'         of freedom, asymptotic p-value.}
#'   \item{\code{RMSEA2}, \code{RMSEA.lower}, \code{RMSEA.upper},
#'         \code{RMSEA.CI.level}, \code{RMSEA.p.close}}{RMSEA estimate,
#'         confidence interval, interval level, close-fit test p-value.}
#'   \item{\code{SRMSR}, \code{RMSR}}{Standardized and raw root mean
#'         square residuals.}
#'   \item{\code{CFI}, \code{TLI}, \code{IFI}}{Incremental fit indices.}
#'   \item{\code{McFadden.R2}, \code{Nagelkerke.R2}, \code{CoxSnell.R2},
#'         \code{AldrichNelson.R2}, \code{VeallZimmermann.R2}}{Pseudo-\eqn{R^2}
#'         measures.}
#'   \item{\code{residuals}, \code{residuals.standardized},
#'         \code{correlation.residuals}}{Residual vectors.}
#'   \item{\code{moments.observed}, \code{moments.predicted}}{First- and
#'         second-order marginal moments.}
#'   \item{\code{q3.mean}, \code{q3.abs.max}, \code{q3.adjusted.abs.max},
#'         \code{q3.p95}}{Yen's Q3 local dependence statistics.}
#'   \item{\code{mean.max.posterior}, \code{median.max.posterior},
#'         \code{min.max.posterior}}{Posterior classification diagnostics.}
#'   \item{\code{null.LogLik}, \code{null.M2.value}, \code{null.df}}{Null
#'         (independence) model diagnostics.}
#'   \item{\code{LRT}, \code{LRT.df}, \code{LRT.p}}{Likelihood ratio test
#'         against the null model.}
#' }
#'
#' @references
#' Akaike, H. (1974). A new look at the statistical model identification.
#' \emph{IEEE Transactions on Automatic Control}, 19, 716--723.
#'
#' Bozdogan, H. (1987). Model selection and Akaike's Information Criterion
#' (AIC): The general theory and its analytical extensions.
#' \emph{Psychometrika}, 52, 345--370. \doi{10.1007/BF02294361}
#'
#' Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing model
#' fit. In K. A. Bollen & J. S. Long (Eds.), \emph{Testing structural
#' equation models} (pp. 136--162). Sage.
#'
#'
#' Christensen, K. B., Makransky, G., & Horton, M. C. (2017). Critical
#' values for Yen's Q3. \emph{Applied Psychological Measurement}, 41,
#' 178--194. \doi{10.1177/0146621616677520}
#'
#' Hannan, E. J., & Quinn, B. G. (1979). The determination of the order of
#' an autoregression. \emph{Journal of the Royal Statistical Society:
#' Series B}, 41, 190--195.
#'
#' Hu, L. T., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in
#' covariance structure analysis. \emph{Structural Equation Modeling}, 6,
#' 1--55. \doi{10.1080/10705519909540118}
#'
#' Hurvich, C. M., & Tsai, C. L. (1989). Regression and time series model
#' selection in small samples. \emph{Biometrika}, 76, 297--307.
#' \doi{10.1093/biomet/76.2.297}
#'
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
#' analysis and determination of sample size for covariance structure
#' modeling. \emph{Psychological Methods}, 1, 130--149.
#'
#' Maydeu-Olivares, A., & Joe, H. (2005). Limited- and full-information
#' estimation and goodness-of-fit testing in 2^n contingency tables: A
#' unified framework. \emph{Journal of the American Statistical
#' Association}, 100, 1009--1020.
#'
#' Maydeu-Olivares, A., & Joe, H. (2006). Limited information
#' goodness-of-fit testing in multidimensional contingency tables.
#' \emph{Psychometrika}, 71, 713--732. \doi{10.1007/s11336-005-1295-9}
#'
#' McFadden, D. (1974). Conditional logit analysis of qualitative choice
#' behavior. In P. Zarembka (Ed.), \emph{Frontiers in econometrics}
#' (pp. 105--142). Academic Press.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals
#' of Statistics}, 6, 461--464. \doi{10.1214/aos/1176344136}
#'
#' Sclove, S. L. (1987). Application of model-selection criteria to some
#' problems in multivariate analysis. \emph{Psychometrika}, 52, 333--343.
#'
#' Yen, W. M. (1984). Effects of local item dependence on the fit and
#' equating performance of the three-parameter logistic model. \emph{Applied
#' Psychological Measurement}, 8, 125--145. \doi{10.1177/014662168400800201}
#'
#' @seealso
#' \code{\link{get.fit.index}} for model-specific wrappers,
#' \code{\link{summary.good.of.fit}} for formatted fit tables,
#' \code{\link{print.good.of.fit}} for the print method.
#'
#' @export
good.of.fit <- function(par.vec, loglik.fun, prob.fun, response, npar, pi,
                        pi.fun = NULL,
                        alpha = 0.05,
                        response.type = c("auto", "binary", "polytomous", "nominal", "ordinary"),
                        response.K = NULL,
                        bivariate.groups = NULL,
                        nominal.groups = NULL,
                        exclusive.groups = NULL,
                        n.boot = 1000L, ...) {

  N <- nrow(response)
  I <- ncol(response)
  n.boot <- as.integer(n.boot[1L])
  if (is.na(n.boot) || n.boot < 0L) n.boot <- 0L

  response <- matrix(as.numeric(response), nrow = N, ncol = I)

  observed.K <- apply(response, 2L, function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_integer_)
    as.integer(max(x) + 1L)
  })
  if (is.null(response.K)) {
    K <- observed.K
  } else {
    K <- as.integer(response.K)
    if (length(K) != I) {
      stop("'response.K' must have one entry per response column.", call. = FALSE)
    }
    if (any(!is.na(observed.K) & observed.K > K)) {
      stop("'response.K' is smaller than an observed response category.", call. = FALSE)
    }
  }
  K[is.na(K)] <- 2L
  response.type <- match.arg(response.type)
  if (response.type == "ordinary") {
    response.type <- "auto"
  }
  if (!is.null(exclusive.groups)) {
    nominal.groups <- exclusive.groups
    response.type <- "nominal"
  }

  # Nominal-data branch:
  # response is already expanded into K-1 binary indicators per nominal item.
  # Indicators from the same original item/block are mutually exclusive; they
  # contribute first-order moments and Xi11 covariances, but no bivariate
  # limited-information margins. TIRT and ordinary binary/polytomous models
  # keep using the default binary/polytomous branches.
  if (!is.null(exclusive.groups)) {
    warning("'exclusive.groups' is deprecated; use response.type = 'nominal' with 'nominal.groups'.",
            call. = FALSE)
  }

  prob_result <- prob.fun(par.vec, ...)
  prob_result <- pmin(pmax(prob_result, 1e-50), 1 - 1e-50)

  get_pi <- function(par.vec, n_support, ...) {
    out <- if (is.null(pi.fun)) pi else pi.fun(par.vec, ...)
    if (is.null(out)) {
      out <- rep(1 / n_support, n_support)
    }
    out <- as.vector(out)
    if (length(out) != n_support) {
      stop("'pi' must have length equal to nrow(prob.fun(par.vec)).", call. = FALSE)
    }
    if (anyNA(out) || any(!is.finite(out)) || any(out < 0) || sum(out) <= 0) {
      stop("'pi' must contain finite non-negative weights with a positive sum.", call. = FALSE)
    }
    out / sum(out)
  }
  pi <- get_pi(par.vec, nrow(prob_result), ...)

  if (response.type == "auto") {
    response.type <- if (ncol(prob_result) == ncol(response)) "binary" else "polytomous"
  }
  if (response.type == "binary" && ncol(prob_result) != ncol(response)) {
    stop("response.type = 'binary' requires prob.fun() to return one probability column per response column.", call. = FALSE)
  }
  if (response.type == "polytomous" && ncol(prob_result) != sum(K)) {
    stop("response.type = 'polytomous' requires prob.fun() to return stacked category probabilities.", call. = FALSE)
  }
  if (response.type == "nominal" && ncol(prob_result) != ncol(response)) {
    stop("response.type = 'nominal' requires prob.fun() to return one probability column per expanded indicator.", call. = FALSE)
  }

  if (response.type == "nominal") {
    nominal.groups <- as.integer(nominal.groups)
    if (length(nominal.groups) != I) {
      stop("'nominal.groups' must have one entry per response column.", call. = FALSE)
    }
    if (any(K > 2L)) {
      stop("response.type = 'nominal' requires expanded binary indicator responses.", call. = FALSE)
    }
    if (is.null(bivariate.groups)) {
      bivariate.groups <- nominal.groups
    }
  }
  if (!is.null(bivariate.groups)) {
    bivariate.groups <- as.integer(bivariate.groups)
    if (length(bivariate.groups) != I) {
      stop("'bivariate.groups' must have one entry per response column.", call. = FALSE)
    }
  }
  is_nominal <- response.type == "nominal"
  is_polytomous <- response.type == "polytomous"

  if (is_polytomous) {
    EI_list <- cpp_compute_EIs(prob_result, K)
    EIs     <- EI_list$EIs
    EIs2    <- EI_list$EIs2
  } else {
    EIs  <- prob_result
    EIs2 <- prob_result
  }

  # Filter out columns with zero observations (e.g. unobserved pairs in MOLE/PICK)
  col.has.obs <- colSums(!is.na(response)) > 0
  prob.col.has.obs <- NULL
  if (!all(col.has.obs)) {
    if (is_polytomous) {
      prob.col.has.obs <- rep(col.has.obs, times = K)
    } else {
      prob.col.has.obs <- col.has.obs
    }
    response <- response[, col.has.obs, drop = FALSE]
    EIs      <- EIs[, col.has.obs, drop = FALSE]
    EIs2     <- EIs2[, col.has.obs, drop = FALSE]
    K        <- K[col.has.obs]
    if (is_nominal) {
      nominal.groups <- nominal.groups[col.has.obs]
    }
    if (!is.null(bivariate.groups)) {
      bivariate.groups <- bivariate.groups[col.has.obs]
    }
    I        <- ncol(response)
  }

  observed.responses.n <- sum(!is.na(response))
  logLik <- loglik.fun(par.vec, ...)

  deviance <- -2 * logLik
  AIC <- deviance + 2 * npar
  AICc <- if (N > (npar + 1L)) {
    AIC + (2 * npar * (npar + 1L)) / (N - npar - 1L)
  } else {
    NA_real_
  }
  BIC <- deviance + log(N) * npar
  CAIC <- deviance + npar * (log(N) + 1)
  SABIC <- deviance + npar * log((N + 2) / 24)
  HQIC <- if (N > exp(1)) deviance + 2 * npar * log(log(N)) else NA_real_
  logLik.per.case <- logLik / N
  logLik.per.response <- if (observed.responses.n > 0L) logLik / observed.responses.n else NA_real_
  neg.logLik.per.case <- -logLik.per.case
  neg.logLik.per.response <- -logLik.per.response
  deviance.per.case <- deviance / N
  deviance.per.response <- if (observed.responses.n > 0L) deviance / observed.responses.n else NA_real_
  perplexity.per.case <- .gof_safe_exp(neg.logLik.per.case)
  perplexity.per.response <- .gof_safe_exp(neg.logLik.per.response)
  RMSEA.CI.level <- 1 - 2 * alpha

  # Compute all bivariate pairs
  pairs.total <- if (I > 1L) cpp_make_pairs(as.integer(I)) else matrix(integer(0), nrow = 2L, ncol = 0L)
  if (is_nominal && I > 1L && !is.null(pairs.total)) {
    pairs.total <- cpp_filter_pairs_between_groups(pairs.total, nominal.groups)
  }
  if (!is_nominal && !is.null(bivariate.groups) && I > 1L && !is.null(pairs.total)) {
    pairs.total <- cpp_filter_pairs_between_groups(pairs.total, bivariate.groups)
  }
  pairs.n <- ncol(pairs.total)

  moments.n <- I + pairs.n
  npar.per.case <- npar / N
  npar.per.moment <- npar / moments.n

  moments.predicted.func <- function(par.vec, ...) {
    prob <- prob.fun(par.vec, ...)
    prob <- pmin(pmax(prob, 1e-50), 1 - 1e-50)
    pi.current <- get_pi(par.vec, nrow(prob), ...)
    if (!is.null(prob.col.has.obs)) {
      prob <- prob[, prob.col.has.obs, drop = FALSE]
    }
    if (is_nominal) {
      return(cpp_moments_predicted_grouped(prob, pi.current, pairs.total, FALSE, K, nominal.groups))
    }
    cpp_moments_predicted(prob, pi.current, pairs.total, is_polytomous, K)
  }

  moments.predicted <- moments.predicted.func(par.vec, ...)

  nq <- nrow(EIs)
  pa_vec <- colSums(EIs * matrix(pi, nq, ncol(EIs), byrow=FALSE))

  # Xi via C++
  Xi_list <- if (is_nominal) {
    cpp_compute_Xi_grouped(EIs, EIs2, pi, pairs.total, nominal.groups)
  } else {
    cpp_compute_Xi(EIs, EIs2, pi, pairs.total)
  }
  Xi2 <- Xi_list$Xi2

  # Jacobian of predicted moments with respect to free parameters
  delta <- numDeriv::jacobian(moments.predicted.func, par.vec, ...)

  observed.summary <- cpp_gof_observed_summary(response, pairs.total)
  moments.observed <- observed.summary$moments
  cross <- observed.summary$cross
  cN    <- observed.summary$count
  mu    <- observed.summary$means
  R.obs <- stats::cov2cor(cross / cN - outer(mu, mu))

  E1 <- pa_vec
  E2 <- if (is_nominal) {
    cpp_compute_E2_grouped(EIs, EIs2, pi, nominal.groups)
  } else {
    cpp_compute_E2(EIs, EIs2, pi)
  }
  Kr <- stats::cov2cor(E2 - outer(E1, E1))
  lt <- lower.tri(R.obs)
  if (!is.null(bivariate.groups)) {
    lt <- lt & outer(bivariate.groups, bivariate.groups, "!=")
  }
  lt.n <- sum(lt, na.rm = TRUE)
  SRMSR <- if (lt.n > 0L) {
    sqrt(sum((R.obs[lt] - Kr[lt])^2, na.rm = TRUE) / lt.n)
  } else {
    NA_real_
  }
  correlation.residuals <- if (lt.n > 0L) R.obs[lt] - Kr[lt] else numeric(0)

  residuals <- moments.observed - moments.predicted
  residual.variance <- diag(Xi2)
  residuals.standardized <- sqrt(N) * residuals / sqrt(residual.variance)
  residuals.standardized[!is.finite(residuals.standardized)] <- NA_real_
  residuals.univariate <- residuals[seq_len(I)]
  residuals.bivariate <- if (pairs.n > 0L) residuals[I + seq_len(pairs.n)] else numeric(0)
  residuals.standardized.univariate <- residuals.standardized[seq_len(I)]
  residuals.standardized.bivariate <- if (pairs.n > 0L) {
    residuals.standardized[I + seq_len(pairs.n)]
  } else {
    numeric(0)
  }

  covariance.residuals <- numeric(0)
  if (pairs.n > 0L) {
    pair.idx1 <- pairs.total[1L, ]
    pair.idx2 <- pairs.total[2L, ]
    obs.pair.moments <- moments.observed[I + seq_len(pairs.n)]
    pred.pair.moments <- moments.predicted[I + seq_len(pairs.n)]
    obs.cov.pairs <- obs.pair.moments -
      moments.observed[pair.idx1] * moments.observed[pair.idx2]
    pred.cov.pairs <- pred.pair.moments - pa_vec[pair.idx1] * pa_vec[pair.idx2]
    covariance.residuals <- obs.cov.pairs - pred.cov.pairs
  }

  finite_summary <- function(x, fun) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    fun(x)
  }
  safe_ratio <- function(numerator, denominator) {
    if (!is.finite(numerator) || !is.finite(denominator) ||
        abs(denominator) < .Machine$double.eps) {
      return(NA_real_)
    }
    numerator / denominator
  }

  RMSR <- finite_summary(residuals, function(x) sqrt(mean(x^2)))
  RMSR.univariate <- finite_summary(residuals.univariate, function(x) sqrt(mean(x^2)))
  RMSR.bivariate <- finite_summary(residuals.bivariate, function(x) sqrt(mean(x^2)))
  mean.abs.residual <- finite_summary(abs(residuals), mean)
  mean.abs.residual.univariate <- finite_summary(abs(residuals.univariate), mean)
  mean.abs.residual.bivariate <- finite_summary(abs(residuals.bivariate), mean)
  max.abs.residual <- finite_summary(abs(residuals), max)
  max.abs.residual.univariate <- finite_summary(abs(residuals.univariate), max)
  max.abs.residual.bivariate <- finite_summary(abs(residuals.bivariate), max)
  RMS.standardized.residual <- finite_summary(residuals.standardized, function(x) sqrt(mean(x^2)))
  RMS.standardized.residual.univariate <- finite_summary(
    residuals.standardized.univariate, function(x) sqrt(mean(x^2))
  )
  RMS.standardized.residual.bivariate <- finite_summary(
    residuals.standardized.bivariate, function(x) sqrt(mean(x^2))
  )
  mean.abs.standardized.residual <- finite_summary(abs(residuals.standardized), mean)
  mean.abs.standardized.residual.univariate <- finite_summary(
    abs(residuals.standardized.univariate), mean
  )
  mean.abs.standardized.residual.bivariate <- finite_summary(
    abs(residuals.standardized.bivariate), mean
  )
  max.abs.standardized.residual <- finite_summary(abs(residuals.standardized), max)
  max.abs.standardized.residual.univariate <- finite_summary(
    abs(residuals.standardized.univariate), max
  )
  max.abs.standardized.residual.bivariate <- finite_summary(
    abs(residuals.standardized.bivariate), max
  )
  prop.abs.standardized.residual.gt.1.96 <- finite_summary(abs(residuals.standardized), function(x) mean(x > 1.96))
  prop.abs.standardized.residual.gt.2.58 <- finite_summary(abs(residuals.standardized), function(x) mean(x > 2.58))
  prop.abs.standardized.residual.gt.3.29 <- finite_summary(abs(residuals.standardized), function(x) mean(x > 3.29))
  mean.abs.correlation.residual <- finite_summary(abs(correlation.residuals), mean)
  max.abs.correlation.residual <- finite_summary(abs(correlation.residuals), max)
  q95.abs.correlation.residual <- finite_summary(abs(correlation.residuals), function(x) {
    as.numeric(stats::quantile(x, probs = 0.95, names = FALSE, type = 8))
  })
  SRMSCR <- finite_summary(covariance.residuals, function(x) sqrt(mean(x^2)))
  mean.abs.covariance.residual <- finite_summary(abs(covariance.residuals), mean)
  max.abs.covariance.residual <- finite_summary(abs(covariance.residuals), max)

  prob.posterior <- prob_result
  if (!is.null(prob.col.has.obs)) {
    prob.posterior <- prob.posterior[, prob.col.has.obs, drop = FALSE]
  }
  posterior.diagnostics <- .gof_posterior_diagnostics(
    prob = prob.posterior,
    response = response,
    pi = pi,
    response.type = response.type,
    K = K,
    nominal.groups = if (is_nominal) nominal.groups else NULL,
    EIs = EIs,
    pair.mask = lt
  )
  null.logLik <- function() {
    cpp_gof_null_loglik(
      response = response,
      K = K,
      response_type = response.type,
      nominal_groups = if (is_nominal) nominal.groups else integer(0),
      first_moments = moments.observed[1:I]
    )
  }

  logLik.nul <- null.logLik()
  deviance.null <- -2 * logLik.nul
  pseudo.R2 <- 1 - safe_ratio(logLik, logLik.nul)
  McFadden.R2 <- pseudo.R2
  McFadden.adj.R2 <- 1 - safe_ratio(logLik - npar, logLik.nul)
  CoxSnell.R2 <- if (is.finite(logLik.nul) && is.finite(logLik)) {
    1 - .gof_safe_exp((2 / N) * (logLik.nul - logLik))
  } else {
    NA_real_
  }
  Nagelkerke.R2 <- if (is.finite(CoxSnell.R2) && is.finite(logLik.nul)) {
    safe_ratio(CoxSnell.R2, 1 - .gof_safe_exp((2 / N) * logLik.nul))
  } else {
    NA_real_
  }

  LRT <- deviance.null - deviance
  LRT.df <- npar
  LRT.p <- if (LRT > 0 && LRT.df > 0L) 1 - stats::pchisq(LRT, LRT.df) else NA_real_
  G2.null <- LRT

  AldrichNelson.R2 <- if (is.finite(LRT) && is.finite(N) && LRT + N > 0)
    LRT / (LRT + N) else NA_real_
  VeallZimmermann.R2 <- if (is.finite(AldrichNelson.R2) && is.finite(logLik.nul))
    AldrichNelson.R2 * ((-2 * logLik.nul + N) / (-2 * logLik.nul)) else NA_real_

  first.mean.null <- moments.observed[1:I]
  first.second.null <- observed.summary$seconds
  first.var.null <- pmax(.Machine$double.xmin, first.second.null - first.mean.null^2)
  diag.Xi2.null <- first.var.null
  if (pairs.n > 0L) {
    pair.mean.null <- first.mean.null[pairs.total[1L, ]] * first.mean.null[pairs.total[2L, ]]
    pair.second.null <- first.second.null[pairs.total[1L, ]] * first.second.null[pairs.total[2L, ]]
    pair.var.null <- pmax(.Machine$double.xmin, pair.second.null - pair.mean.null^2)
    diag.Xi2.null <- c(diag.Xi2.null, pair.var.null)
  }

  moments.observed.null <- numeric(moments.n)
  moments.observed.null[1:I] <- moments.observed[1:I]
  if (I > 1L) {
    for (k in seq_len(pairs.n)) {
      ii <- pairs.total[1L, k]
      jj <- pairs.total[2L, k]
      moments.observed.null[I + k] <- moments.observed.null[ii] * moments.observed.null[jj]
    }
  }

  residuals.null <- moments.observed - moments.observed.null
  M2.null <- as.numeric(N * sum(residuals.null^2 / diag.Xi2.null))
  df.null <- moments.n

  extended.results <- list(
    logLik.per.response = logLik.per.response,
    neg.logLik.per.case = neg.logLik.per.case,
    neg.logLik.per.response = neg.logLik.per.response,
    deviance.per.case = deviance.per.case,
    deviance.per.response = deviance.per.response,
    perplexity.per.case = perplexity.per.case,
    perplexity.per.response = perplexity.per.response,
    observed.responses.n = observed.responses.n,
    npar.per.case = npar.per.case,
    npar.per.moment = npar.per.moment,
    RMSR.univariate = RMSR.univariate,
    RMSR.bivariate = RMSR.bivariate,
    mean.abs.residual.univariate = mean.abs.residual.univariate,
    mean.abs.residual.bivariate = mean.abs.residual.bivariate,
    max.abs.residual.univariate = max.abs.residual.univariate,
    max.abs.residual.bivariate = max.abs.residual.bivariate,
    RMS.standardized.residual.univariate = RMS.standardized.residual.univariate,
    RMS.standardized.residual.bivariate = RMS.standardized.residual.bivariate,
    mean.abs.standardized.residual.univariate = mean.abs.standardized.residual.univariate,
    mean.abs.standardized.residual.bivariate = mean.abs.standardized.residual.bivariate,
    max.abs.standardized.residual.univariate = max.abs.standardized.residual.univariate,
    max.abs.standardized.residual.bivariate = max.abs.standardized.residual.bivariate,
    q95.abs.correlation.residual = q95.abs.correlation.residual,
    SRMSCR = SRMSCR,
    mean.abs.covariance.residual = mean.abs.covariance.residual,
    max.abs.covariance.residual = max.abs.covariance.residual,
    pseudo.R2 = pseudo.R2,
    McFadden.R2 = McFadden.R2,
    McFadden.adj.R2 = McFadden.adj.R2,
    CoxSnell.R2 = CoxSnell.R2,
    Nagelkerke.R2 = Nagelkerke.R2,
    AldrichNelson.R2 = AldrichNelson.R2,
    VeallZimmermann.R2 = VeallZimmermann.R2,
    LRT = LRT,
    LRT.df = LRT.df,
    LRT.p = LRT.p,
    G2.null = G2.null,
    case.logLik = posterior.diagnostics$case.logLik,
    case.logLik.summary = posterior.diagnostics$case.logLik.summary,
    mean.max.posterior = posterior.diagnostics$mean.max.posterior,
    median.max.posterior = posterior.diagnostics$median.max.posterior,
    min.max.posterior = posterior.diagnostics$min.max.posterior,
    q3.mean = posterior.diagnostics$q3.mean,
    q3.max = posterior.diagnostics$q3.max,
    q3.abs.max = posterior.diagnostics$q3.abs.max,
    q3.p95 = posterior.diagnostics$q3.p95,
    q3.adjusted.max = posterior.diagnostics$q3.adjusted.max,
    q3.adjusted.abs.max = posterior.diagnostics$q3.adjusted.abs.max,
    q3.adjusted.p95 = posterior.diagnostics$q3.adjusted.p95,
    covariance.residuals = covariance.residuals
  )

  delta.qr <- qr(delta)
  if ((ncol(delta) + 1L) > ncol(qr.Q(delta.qr, complete = TRUE))) {
    warning('M2 cannot be calculated: too few degrees of freedom. M2 and related indices set to NA.', call. = FALSE)

    results <- list(
      LogLik            = logLik,
      logLik.per.case   = logLik.per.case,
      deviance          = deviance,
      npar              = npar,
      N                 = N,
      nitems            = I,
      response.type     = response.type,
      moments.n         = moments.n,
      pairs.n           = pairs.n,
      AIC               = AIC,
      AICc              = AICc,
      BIC               = BIC,
      CAIC              = CAIC,
      SABIC             = SABIC,
      HQIC              = HQIC,
      M2.value          = NA_real_,
      df                = NA_real_,
      p.value           = NA_real_,
      RMSEA2            = NA_real_,
      RMSEA.lower       = NA_real_,
      RMSEA.upper       = NA_real_,
      RMSEA.CI.level    = RMSEA.CI.level,
      RMSEA.p.close     = NA_real_,
      SRMSR             = SRMSR,
      RMSR              = RMSR,
      mean.abs.residual = mean.abs.residual,
      max.abs.residual  = max.abs.residual,
      RMS.standardized.residual = RMS.standardized.residual,
      mean.abs.standardized.residual = mean.abs.standardized.residual,
      max.abs.standardized.residual = max.abs.standardized.residual,
      prop.abs.standardized.residual.gt.1.96 = prop.abs.standardized.residual.gt.1.96,
      prop.abs.standardized.residual.gt.2.58 = prop.abs.standardized.residual.gt.2.58,
      prop.abs.standardized.residual.gt.3.29 = prop.abs.standardized.residual.gt.3.29,
      mean.abs.correlation.residual = mean.abs.correlation.residual,
      max.abs.correlation.residual = max.abs.correlation.residual,
      pseudo.R2               = pseudo.R2,
      AldrichNelson.R2        = AldrichNelson.R2,
      VeallZimmermann.R2      = VeallZimmermann.R2,
      CFI               = NA_real_,
      TLI               = NA_real_,
      IFI               = NA_real_,
      McDonald.NCI      = NA_real_,
      null.LogLik       = logLik.nul,
      null.deviance     = deviance.null,
      null.M2.value     = M2.null,
      null.df           = df.null,
      LRT               = LRT,
      LRT.df            = LRT.df,
      LRT.p             = LRT.p,
      G2.null           = G2.null,
      residuals         = residuals,
      residuals.standardized = residuals.standardized,
      correlation.residuals = correlation.residuals,
      moments.observed  = moments.observed,
      moments.predicted = moments.predicted
    )

    for (nm in names(extended.results)) {
      results[[nm]] <- extended.results[[nm]]
    }
    class(results) <- "good.of.fit"
    return(results)
  }
  deltac <- qr.Q(delta.qr, complete = TRUE)[, (ncol(delta) + 1L):ncol(qr.Q(delta.qr, complete = TRUE)), drop = FALSE]

  inner <- crossprod(deltac, Xi2 %*% deltac)
  C2 <- tryCatch(
    deltac %*% solve(inner) %*% t(deltac),
    error = function(e) {
      deltac %*% MASS::ginv(inner) %*% t(deltac)
    }
  )

  M2.value <- as.numeric(abs(N * crossprod(residuals, C2 %*% residuals)))

  df <- moments.n - npar
  p.value <- 1 - stats::pchisq(M2.value, df)
  RMSEA2 <- sqrt(max(0, M2.value - df) / (N * df))
  RMSEA2.CI.value <- RMSEA.CI(M2.value, df, N, ci.lower = alpha, ci.upper = 1-alpha)
  RMSEA.lower <- RMSEA2.CI.value[1]
  RMSEA.upper <- RMSEA2.CI.value[2]
  RMSEA.p.close <- 1 - stats::pchisq(
    M2.value, df = df, ncp = N * df * 0.05^2
  )
  McDonald.NCI <- exp(-max(M2.value - df, 0) / (2 * N))

  CFI <- NA_real_; TLI <- NA_real_
  if (M2.null > M2.value) {
    CFI <- min(1, max(0, 1 - (M2.value - df) / (M2.null - df.null)))
    TLI <- (M2.null / df.null - M2.value / df) / (M2.null / df.null - 1)
  } else {
    warning('Null model M2 <= fitted model M2; CFI/TLI not computed.', call. = FALSE)
  }
  IFI <- safe_ratio(M2.null - M2.value, M2.null - df)

  results <- list(
    LogLik            = logLik,
    logLik.per.case   = logLik.per.case,
    deviance          = deviance,
    npar              = npar,
    N                 = N,
    nitems            = I,
    response.type     = response.type,
    moments.n         = moments.n,
    pairs.n           = pairs.n,
    AIC               = AIC,
    AICc              = AICc,
    BIC               = BIC,
    CAIC              = CAIC,
    SABIC             = SABIC,
    HQIC              = HQIC,
    M2.value          = M2.value,
    df                = df,
    p.value           = p.value,
    RMSEA2            = RMSEA2,
    RMSEA.lower       = RMSEA.lower,
    RMSEA.upper       = RMSEA.upper,
    RMSEA.CI.level    = RMSEA.CI.level,
    RMSEA.p.close     = RMSEA.p.close,
    SRMSR             = SRMSR,
    RMSR              = RMSR,
    mean.abs.residual = mean.abs.residual,
    max.abs.residual  = max.abs.residual,
    RMS.standardized.residual = RMS.standardized.residual,
    mean.abs.standardized.residual = mean.abs.standardized.residual,
    max.abs.standardized.residual = max.abs.standardized.residual,
    prop.abs.standardized.residual.gt.1.96 = prop.abs.standardized.residual.gt.1.96,
    prop.abs.standardized.residual.gt.2.58 = prop.abs.standardized.residual.gt.2.58,
    prop.abs.standardized.residual.gt.3.29 = prop.abs.standardized.residual.gt.3.29,
    mean.abs.correlation.residual = mean.abs.correlation.residual,
    max.abs.correlation.residual = max.abs.correlation.residual,
    pseudo.R2               = pseudo.R2,
    AldrichNelson.R2        = AldrichNelson.R2,
    VeallZimmermann.R2      = VeallZimmermann.R2,
    CFI               = CFI,
    TLI               = TLI,
    IFI               = IFI,
    McDonald.NCI      = McDonald.NCI,
    null.LogLik       = logLik.nul,
    null.deviance     = deviance.null,
    null.M2.value     = M2.null,
    null.df           = df.null,
    LRT               = LRT,
    LRT.df            = LRT.df,
    LRT.p             = LRT.p,
    G2.null           = G2.null,
    residuals         = residuals,
    residuals.standardized = residuals.standardized,
    correlation.residuals = correlation.residuals,
    moments.observed  = moments.observed,
    moments.predicted = moments.predicted
  )

  for (nm in names(extended.results)) {
    results[[nm]] <- extended.results[[nm]]
  }
  class(results) <- "good.of.fit"

  return(results)
}

RMSEA.CI <- function(M2, df, N, ci.lower = 0.05, ci.upper = 0.95){

  lower.lambda <- function(lambda) {
    stats::pchisq(M2, df = df, ncp = lambda) - ci.upper
  }
  upper.lambda <- function(lambda) {
    stats::pchisq(M2, df = df, ncp = lambda) - ci.lower
  }

  lambda.l <- try(
    stats::uniroot(f = lower.lambda, lower = 0, upper = M2)$root,
    silent = TRUE
  )
  lambda.u <- try(
    stats::uniroot(
      f = upper.lambda, lower = 0, upper = max(N, M2 * 5)
    )$root,
    silent = TRUE
  )

  if (!is(lambda.l, "try-error")) {
    RMSEA.lower <- sqrt(lambda.l/(N * df))
  }else {
    RMSEA.lower <- 0
  }

  if (!is(lambda.u, "try-error")) {
    RMSEA.upper <- sqrt(lambda.u/(N * df))
  }else {
    RMSEA.upper <- 0
  }

  return(c(RMSEA.lower, RMSEA.upper))
}

.gof_safe_exp <- function(x) {
  out <- rep(NA_real_, length(x))
  finite <- is.finite(x)
  out[finite] <- exp(pmin(x[finite], log(.Machine$double.xmax)))
  out[is.infinite(x) & x < 0] <- 0
  if (length(out) == 1L) out[[1L]] else out
}

.gof_log_sum_exp <- function(x) {
  m <- max(x, na.rm = TRUE)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

.gof_summary6 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(Min = NA_real_, Q1 = NA_real_, Median = NA_real_,
             Mean = NA_real_, Q3 = NA_real_, Max = NA_real_))
  }
  c(
    Min = min(x),
    Q1 = as.numeric(stats::quantile(x, 0.25, names = FALSE, type = 8)),
    Median = stats::median(x),
    Mean = mean(x),
    Q3 = as.numeric(stats::quantile(x, 0.75, names = FALSE, type = 8)),
    Max = max(x)
  )
}

.gof_empty_posterior_diagnostics <- function() {
  list(
    case.logLik = numeric(0),
    case.logLik.summary = .gof_summary6(numeric(0)),
    mean.max.posterior = NA_real_,
    median.max.posterior = NA_real_,
    min.max.posterior = NA_real_,
    q3.mean = NA_real_,
    q3.max = NA_real_,
    q3.abs.max = NA_real_,
    q3.p95 = NA_real_,
    q3.adjusted.max = NA_real_,
    q3.adjusted.abs.max = NA_real_,
    q3.adjusted.p95 = NA_real_
  )
}

.gof_case_loglik <- function(prob, response, pi, response.type, K,
                             nominal.groups = NULL) {
  eps <- 1e-50
  prob <- pmin(pmax(prob, eps), 1 - eps)
  pi <- as.vector(pi)
  pi <- pmax(pi, eps)
  pi <- pi / sum(pi)
  N <- nrow(response)
  nq <- nrow(prob)
  I <- ncol(response)
  log.joint <- matrix(log(pi), nrow = N, ncol = nq, byrow = TRUE)

  if (response.type == "binary") {
    log.p <- log(prob)
    log.q <- log1p(-prob)
    for (j in seq_len(I)) {
      y <- response[, j]
      idx1 <- which(!is.na(y) & y >= 0.5)
      idx0 <- which(!is.na(y) & y < 0.5)
      if (length(idx1)) {
        log.joint[idx1, ] <- log.joint[idx1, , drop = FALSE] +
          matrix(log.p[, j], nrow = length(idx1), ncol = nq, byrow = TRUE)
      }
      if (length(idx0)) {
        log.joint[idx0, ] <- log.joint[idx0, , drop = FALSE] +
          matrix(log.q[, j], nrow = length(idx0), ncol = nq, byrow = TRUE)
      }
    }
  } else if (response.type == "polytomous") {
    offsets <- cumsum(c(0L, K[-length(K)]))
    for (j in seq_len(I)) {
      y <- response[, j]
      bad <- which(!is.na(y) & (y < 0 | y >= K[j] | y != floor(y)))
      if (length(bad)) {
        log.joint[bad, ] <- -Inf
      }
      for (cat in seq_len(K[j]) - 1L) {
        idx <- which(!is.na(y) & y == cat)
        if (length(idx)) {
          col <- offsets[j] + cat + 1L
          log.joint[idx, ] <- log.joint[idx, , drop = FALSE] +
            matrix(log(prob[, col]), nrow = length(idx), ncol = nq, byrow = TRUE)
        }
      }
    }
  } else if (response.type == "nominal") {
    groups <- unique(nominal.groups)
    for (g in groups) {
      cols <- which(nominal.groups == g)
      yy <- response[, cols, drop = FALSE]
      complete <- rowSums(is.na(yy)) == 0L
      active <- rowSums(yy == 1, na.rm = TRUE)
      invalid <- which(complete & active > 1L)
      if (length(invalid)) {
        log.joint[invalid, ] <- -Inf
      }
      idx0 <- which(complete & active == 0L)
      if (length(idx0)) {
        p0 <- pmax(eps, 1 - rowSums(prob[, cols, drop = FALSE]))
        log.joint[idx0, ] <- log.joint[idx0, , drop = FALSE] +
          matrix(log(p0), nrow = length(idx0), ncol = nq, byrow = TRUE)
      }
      for (h in seq_along(cols)) {
        idx <- which(complete & active == 1L & yy[, h] == 1)
        if (length(idx)) {
          log.joint[idx, ] <- log.joint[idx, , drop = FALSE] +
            matrix(log(prob[, cols[h]]), nrow = length(idx), ncol = nq, byrow = TRUE)
        }
      }
    }
  } else {
    stop("Unsupported response.type for posterior diagnostics.", call. = FALSE)
  }

  case.logLik <- apply(log.joint, 1L, .gof_log_sum_exp)
  posterior <- matrix(NA_real_, nrow = N, ncol = nq)
  valid <- is.finite(case.logLik)
  if (any(valid)) {
    posterior[valid, ] <- exp(sweep(log.joint[valid, , drop = FALSE],
                                    1L, case.logLik[valid], "-"))
  }
  list(case.logLik = case.logLik, posterior = posterior)
}

.gof_posterior_diagnostics <- function(prob, response, pi, response.type, K,
                                       nominal.groups, EIs, pair.mask) {
  out <- .gof_empty_posterior_diagnostics()
  case <- tryCatch(
    .gof_case_loglik(prob, response, pi, response.type, K, nominal.groups),
    error = function(e) NULL
  )
  if (is.null(case)) {
    return(out)
  }

  posterior <- case$posterior
  valid <- rowSums(is.na(posterior)) == 0L
  out$case.logLik <- case$case.logLik
  out$case.logLik.summary <- .gof_summary6(case$case.logLik)
  if (!any(valid)) {
    return(out)
  }

  post.valid <- posterior[valid, , drop = FALSE]
  max.post <- apply(post.valid, 1L, max)
  out$mean.max.posterior <- mean(max.post)
  out$median.max.posterior <- stats::median(max.post)
  out$min.max.posterior <- min(max.post)

  if (ncol(response) > 1L && ncol(EIs) == ncol(response)) {
    expected <- posterior %*% EIs
    residual.matrix <- response - expected
    residual.matrix[is.na(response)] <- NA_real_
    q3 <- suppressWarnings(stats::cor(residual.matrix, use = "pairwise.complete.obs"))
    q3.values <- q3[pair.mask]
    q3.values <- q3.values[is.finite(q3.values)]
    if (length(q3.values)) {
      q3.adjusted <- q3.values - mean(q3.values)
      out$q3.mean <- mean(q3.values)
      out$q3.max <- max(q3.values)
      out$q3.abs.max <- max(abs(q3.values))
      out$q3.p95 <- as.numeric(stats::quantile(q3.values, 0.95, names = FALSE, type = 8))
      out$q3.adjusted.max <- max(q3.adjusted)
      out$q3.adjusted.abs.max <- max(abs(q3.adjusted))
      out$q3.adjusted.p95 <- as.numeric(stats::quantile(q3.adjusted, 0.95, names = FALSE, type = 8))
    }
  }

  out
}
