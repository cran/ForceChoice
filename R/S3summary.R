#' @title S3 Methods: summary
#'
#' @description
#' Generates structured, comprehensive summaries of objects produced by the
#' \pkg{ForceChoice} package. This generic function dispatches to
#' class-specific methods that extract and organize key information including
#' model configurations, fit statistics, parameter estimates, and convergence
#' diagnostics. Designed for programmatic access and downstream reporting.
#'
#' @param object An object of one of the following classes:
#'   \itemize{
#'     \item Fitted model objects: \code{"MIRT"}, \code{"MGPCM"},
#'       \code{"MGGUM"}, \code{"FCMIRT"}, \code{"FCDCM"}, \code{"FCGGUM"},
#'       \code{"FCGDINA"}, \code{"TIRT"}
#'     \item Goodness-of-fit objects: \code{"good.of.fit"}
#'   }
#' @param digits Number of decimal places for numeric output (default: 4).
#'   Applied uniformly across all methods.
#' @param ... Additional arguments passed to or from other methods (currently
#'   ignored).
#'
#' @return Invisibly returns a structured list containing summary components.
#'   The exact structure depends on the class of \code{object}. All returned
#'   objects carry an appropriate S3 class (e.g., \code{"summary.MIRT"},
#'   \code{"summary.good.of.fit"}) for use with corresponding \code{print}
#'   methods.
#'
#' @details
#' Each method returns a named list with class-specific components optimized
#' for structured access and formatted printing:
#'
#' \describe{
#'   \item{\strong{Fitted Model Objects}}{
#'     Returns a \code{summary.<Class>} object with components:
#'     \describe{
#'       \item{\code{call}}{Original function call.}
#'       \item{\code{model.info}}{List: \code{family}, \code{description},
#'         \code{D}, \code{method}, and model-specific fields (\code{model},
#'         \code{fc.type}, \code{dcm.type}, \code{I.states}, \code{N.block}).}
#'       \item{\code{data.info}}{List: \code{N}, \code{I} or \code{B},
#'         \code{response.type}.}
#'       \item{\code{fit.stats}}{List: \code{LogLik}, \code{npar},
#'         \code{AIC}, \code{BIC}.}
#'       \item{\code{par.summary}}{Matrix: parameter summary (Min, Q1, Median,
#'         Mean, Q3, Max) across items/statements/blocks for each parameter type.}
#'       \item{\code{theta.summary}}{Matrix: person parameter summary per
#'         dimension.}
#'       \item{\code{Corr}}{Matrix: factor inter-trait correlation matrix
#'         (where applicable).}
#'       \item{\code{convergence}}{List: algorithm, batch counts, convergence
#'         flags, latent grid length \code{L} (iStEM) or diagnostic note (Stan).}
#'       \item{\code{digits}}{Numeric: precision used for formatting.}
#'     }
#'     Additional model-specific components:
#'     \describe{
#'       \item{\code{delta.summary}}{FCDCM only: higher-order delta parameters.}
#'       \item{\code{classification}}{FCDCM only: observed attribute patterns,
#'         mean maximum posterior.}
#'       \item{\code{lambda.summary}}{TIRT only: factor loading summary.}
#'     }
#'   }
#'
#'   \item{\strong{Goodness-of-Fit Objects (\code{good.of.fit})}}{
#'     Returns a \code{summary.good.of.fit} object with components:
#'     \describe{
#'       \item{\code{Sample_Size}, \code{Items}, \code{Response_Type},
#'         \code{Parameters}, \code{Moments}, \code{LogLik},
#'         \code{Deviance}}{Overview statistics.}
#'       \item{\code{IC}}{Named vector: AIC, AICc, BIC, CAIC, SABIC, HQIC.}
#'       \item{\code{M2}}{List: M2 statistic, df, p-value, RMSEA with CI.}
#'       \item{\code{Residual_Fit}}{Named vector: SRMSR, RMSR, max
#'         standardized residual.}
#'       \item{\code{Comparative_Fit}}{Named vector: CFI, TLI, IFI.}
#'       \item{\code{Pseudo_R2}}{Named vector: McFadden, Cox--Snell,
#'         Nagelkerke, Aldrich--Nelson, Veall--Zimmermann.}
#'       \item{\code{Local_Dependence}}{Named vector: Q3 mean, max, adjusted
#'         max, P95.}
#'       \item{\code{Classification}}{Named vector: posterior entropy,
#'         max posterior.}
#'       \item{\code{Null_Model}}{Named vector: null model LogLik, deviance,
#'         M2, df.}
#'       \item{\code{tables}}{List of data frames for formatted printing.}
#'     }
#'   }
#' }
#'
#' @name summary
NULL


# ===========================================================================
# summary.good.of.fit
# ===========================================================================

#' @describeIn summary Summary method for \code{good.of.fit} objects.
#'   Extracts and structures all fit indices: information criteria,
#'   limited-information M2, RMSEA, SRMSR, CFI/TLI/IFI, pseudo-\eqn{R^2},
#'   local dependence (Yen's Q3), posterior classification, and null-model
#'   diagnostics.
#' @method summary good.of.fit
#' @export
#' @rawNamespace export(summary.good.of.fit)
summary.good.of.fit <- function(object, digits = 4, ...) {

  safe_extract <- function(obj, name) {
    if (!is.null(obj[[name]])) obj[[name]] else NA_real_
  }

  res <- list(
    Sample_Size = object$N,
    Items       = safe_extract(object, "nitems"),
    Response_Type = if (!is.null(object$response.type)) object$response.type else NA_character_,
    Parameters  = object$npar,
    Moments     = safe_extract(object, "moments.n"),
    Bivariate_Moments = safe_extract(object, "pairs.n"),
    Observed_Responses = safe_extract(object, "observed.responses.n"),
    LogLik      = object$LogLik,
    LogLik_per_case = safe_extract(object, "logLik.per.case"),
    LogLik_per_response = safe_extract(object, "logLik.per.response"),
    Deviance    = safe_extract(object, "deviance"),

    IC = c(
      AIC   = safe_extract(object, "AIC"),
      AICc  = safe_extract(object, "AICc"),
      BIC   = safe_extract(object, "BIC"),
      CAIC  = safe_extract(object, "CAIC"),
      SABIC = safe_extract(object, "SABIC"),
      HQIC  = safe_extract(object, "HQIC")
    ),

    M2 = list(
      Statistic   = safe_extract(object, "M2.value"),
      df          = safe_extract(object, "df"),
      p_value     = safe_extract(object, "p.value"),
      RMSEA       = safe_extract(object, "RMSEA2"),
      RMSEA_lower = safe_extract(object, "RMSEA.lower"),
      RMSEA_upper = safe_extract(object, "RMSEA.upper"),
      RMSEA_CI_level = safe_extract(object, "RMSEA.CI.level"),
      RMSEA_p_close = safe_extract(object, "RMSEA.p.close")
    ),

    Residual_Fit = c(
      SRMSR = safe_extract(object, "SRMSR"),
      RMSR = safe_extract(object, "RMSR"),
      MaxAbsStdResidual = safe_extract(object, "max.abs.standardized.residual"),
      Q95AbsCorrResidual = safe_extract(object, "q95.abs.correlation.residual")
    ),

    Comparative_Fit = c(
      CFI       = safe_extract(object, "CFI"),
      TLI       = safe_extract(object, "TLI"),
      IFI       = safe_extract(object, "IFI")
    ),

    Pseudo_R2 = c(
      McFadden_R2        = safe_extract(object, "McFadden.R2"),
      McFadden_adj_R2    = safe_extract(object, "McFadden.adj.R2"),
      CoxSnell_R2        = safe_extract(object, "CoxSnell.R2"),
      Nagelkerke_R2      = safe_extract(object, "Nagelkerke.R2"),
      AldrichNelson_R2   = safe_extract(object, "AldrichNelson.R2"),
      VeallZimmermann_R2 = safe_extract(object, "VeallZimmermann.R2")
    ),

    Local_Dependence = c(
      Q3_mean = safe_extract(object, "q3.mean"),
      Q3_max = safe_extract(object, "q3.max"),
      Q3_abs_max = safe_extract(object, "q3.abs.max"),
      Q3_p95 = safe_extract(object, "q3.p95"),
      Adjusted_Q3_max = safe_extract(object, "q3.adjusted.max"),
      Adjusted_Q3_abs_max = safe_extract(object, "q3.adjusted.abs.max"),
      Adjusted_Q3_p95 = safe_extract(object, "q3.adjusted.p95")
    ),

    Classification = c(
      Posterior_entropy = safe_extract(object, "posterior.entropy"),
      Normalized_entropy = safe_extract(object, "posterior.entropy.normalized"),
      Mean_max_posterior = safe_extract(object, "mean.max.posterior"),
      Median_max_posterior = safe_extract(object, "median.max.posterior"),
      Min_max_posterior = safe_extract(object, "min.max.posterior")
    ),

    Null_Model = c(
      LogLik = safe_extract(object, "null.LogLik"),
      Deviance = safe_extract(object, "null.deviance"),
      M2 = safe_extract(object, "null.M2.value"),
      df = safe_extract(object, "null.df")
    )
  )

  res$tables <- .gof_summary_tables(object, digits = digits)

  class(res) <- "summary.good.of.fit"
  attr(res, "digits") <- digits
  return(res)
}


# ===========================================================================
# Internal helpers for model summaries
# ===========================================================================

.fc_summary_theta <- function(theta_est, digits) {
  if (is.null(theta_est)) return(NULL)
  est <- as.matrix(theta_est)
  D <- ncol(est)
  stats <- matrix(NA_real_, nrow = D, ncol = 6)
  rownames(stats) <- paste0("Dim.", seq_len(D))
  colnames(stats) <- c("Min", "Q1", "Median", "Mean", "Q3", "Max")
  for (d in seq_len(D)) {
    x <- est[, d]
    x <- x[is.finite(x)]
    if (length(x)) {
      stats[d, ] <- c(min(x), stats::quantile(x, 0.25), stats::median(x),
                      mean(x), stats::quantile(x, 0.75), max(x))
    }
  }
  round(stats, digits)
}

.fc_summary_par_range <- function(par_est, par_free, digits) {
  if (is.null(par_est)) return(NULL)
  est <- as.matrix(par_est)
  free <- if (is.null(par_free)) matrix(TRUE, nrow(est), ncol(est)) else as.matrix(par_free)
  nc <- ncol(est)
  res <- matrix(NA_real_, nrow = nc, ncol = 6)
  rownames(res) <- colnames(est)
  colnames(res) <- c("Min", "Q1", "Median", "Mean", "Q3", "Max")
  for (j in seq_len(nc)) {
    vals <- est[, j]
    if (!is.null(free) && ncol(free) == nc) {
      vals <- vals[free[, j]]
    }
    vals <- vals[is.finite(vals)]
    if (length(vals) && !all(vals == 0)) {
      res[j, ] <- c(min(vals), stats::quantile(vals, 0.25),
                    stats::median(vals), mean(vals),
                    stats::quantile(vals, 0.75), max(vals))
    }
  }
  keep <- rowSums(is.finite(res)) > 0L
  if (!any(keep)) return(NULL)
  round(res[keep, , drop = FALSE], digits)
}

.fc_summary_corr <- function(Corr_est, digits) {
  if (is.null(Corr_est)) return(NULL)
  corr <- as.matrix(Corr_est)
  round(corr, digits)
}

.fc_summary_convergence <- function(object) {
  if (!is.null(object$iStEM) && is.list(object$iStEM) && length(object$iStEM)) {
    istem <- object$iStEM
    conv <- list(
      algorithm = "iStEM (Improved Stochastic EM)",
      burn_in_batches = istem$burn.in.size,
      total_batches   = istem$total.number.of.batch,
      burnin_converged = istem$burnin.converged,
      final_converged  = if (!is.null(istem$converged)) istem$converged else NA,
      L = if (!is.null(istem$L)) {
        istem$L
      } else {
        NA
      }
    )
    return(conv)
  }
  if (!is.null(object$EM) && is.list(object$EM) && length(object$EM)) {
    em <- object$EM
    boot <- em$bootstrap
    conv <- list(
      algorithm = "EM (posterior-weight deterministic EM)",
      iterations = em$iter,
      final_converged = if (!is.null(em$converged)) em$converged else NA,
      final_logLik = em$final.logLik,
      final_ll_change = em$final.ll.change,
      final_max_delta_change = em$final.max.delta.change,
      bootstrap_R = if (!is.null(boot)) boot$R else 0L,
      bootstrap_successful = if (!is.null(boot)) boot$successful else 0L
    )
    return(conv)
  }
  if (!is.null(object$stan.obj)) {
    conv <- list(
      algorithm = "Stan (Hamiltonian Monte Carlo)",
      note = "See stan.obj for full convergence diagnostics."
    )
    return(conv)
  }
  list(algorithm = object$method, note = "No detailed convergence info available.")
}

.fc_null_default <- function(x, default) {
  if (is.null(x)) default else x
}

.fc_model_name <- function(class_name, object) {
  switch(class_name,
    MIRT   = paste0("Multidimensional IRT (",
                    toupper(.fc_null_default(.fc_null_default(object$arguments$model, object$model), "m2pl")),
                    ")"),
    MGPCM  = "Multidimensional Generalized Partial Credit Model",
    MGGUM  = "Multidimensional Generalized Graded Unfolding Model",
    FCMIRT = "Forced-Choice Multidimensional IRT",
    FCDCM  = "Forced-Choice Diagnostic Classification Model",
    FCGDINA = "Forced-Choice GDINA",
    FCGGUM = "Forced-Choice Generalized Graded Unfolding Model",
    TIRT   = "Thurstonian IRT for Forced-Choice",
    class_name
  )
}


# ===========================================================================
# Model summary methods
# ===========================================================================

# ---- summary.MIRT ----

#' @describeIn summary Summary method for \code{MIRT} objects.
#'   Multidimensional IRT (1PL--4PL): extracts binary-response model
#'   configuration, fit statistics (LogLik, AIC, BIC), factor correlation
#'   matrix, person and item parameter summaries, and convergence
#'   diagnostics.
#' @method summary MIRT
#' @export
summary.MIRT <- function(object, digits = 4, ...) {
  response <- object$arguments$response
  N <- nrow(response)
  I <- ncol(response)
  D <- object$arguments$D
  model <- if (!is.null(object$arguments$model)) object$arguments$model
           else if (!is.null(object$model)) object$model else "m2pl"
  logLik.val <- as.numeric(object$logLik)

  AIC <- -2 * logLik.val + 2 * object$npar
  BIC <- -2 * logLik.val + object$npar * log(N)

  res <- list(
    call       = object$call,
    model.info = list(
      family = "MIRT",
      description = .fc_model_name("MIRT", object),
      model  = model,
      D      = D,
      method = object$method
    ),
    data.info = list(
      N = N,
      I = I,
      response.type = "binary (0/1)"
    ),
    fit.stats = list(
      LogLik = logLik.val,
      npar   = object$npar,
      AIC    = AIC,
      BIC    = BIC
    ),
    par.summary   = .fc_summary_par_range(object$par$est, object$par$free, digits),
    theta.summary = .fc_summary_theta(object$theta$est, digits),
    Corr          = .fc_summary_corr(object$Corr$est, digits),
    convergence   = .fc_summary_convergence(object),
    digits        = digits
  )
  class(res) <- "summary.MIRT"
  invisible(res)
}

# ---- summary.MGPCM ----

#' @describeIn summary Summary method for \code{MGPCM} objects.
#'   Multidimensional Generalized Partial Credit Model: polytomous
#'   responses with category-specific step parameters.
#' @method summary MGPCM
#' @export
summary.MGPCM <- function(object, digits = 4, ...) {
  response <- object$arguments$response
  N <- nrow(response)
  I <- ncol(response)
  D <- object$arguments$D
  logLik.val <- as.numeric(object$logLik)

  AIC <- -2 * logLik.val + 2 * object$npar
  BIC <- -2 * logLik.val + object$npar * log(N)

  cat_range <- if (!is.null(object$length.poly)) {
    paste0(min(object$length.poly), "--", max(object$length.poly))
  } else "unknown"

  res <- list(
    call       = object$call,
    model.info = list(
      family = "MGPCM",
      description = .fc_model_name("MGPCM", object),
      D      = D,
      method = object$method
    ),
    data.info = list(
      N = N,
      I = I,
      response.type = paste0("polytomous (", cat_range, " categories)")
    ),
    fit.stats = list(
      LogLik = logLik.val,
      npar   = object$npar,
      AIC    = AIC,
      BIC    = BIC
    ),
    par.summary   = .fc_summary_par_range(object$par$est, object$par$free, digits),
    theta.summary = .fc_summary_theta(object$theta$est, digits),
    Corr          = .fc_summary_corr(object$Corr$est, digits),
    convergence   = .fc_summary_convergence(object),
    digits        = digits
  )
  class(res) <- "summary.MGPCM"
  invisible(res)
}

# ---- summary.MGGUM ----

#' @describeIn summary Summary method for \code{MGGUM} objects.
#'   Multidimensional Generalized Graded Unfolding Model: ideal-point
#'   polytomous responses with discrimination (a), location (delta), and
#'   threshold (tau) parameters.
#' @method summary MGGUM
#' @export
summary.MGGUM <- function(object, digits = 4, ...) {
  response <- object$arguments$response
  N <- nrow(response)
  I <- ncol(response)
  D <- object$arguments$D
  logLik.val <- as.numeric(object$logLik)

  AIC <- -2 * logLik.val + 2 * object$npar
  BIC <- -2 * logLik.val + object$npar * log(N)

  cat_range <- if (!is.null(object$length.poly)) {
    paste0(min(object$length.poly), "--", max(object$length.poly))
  } else "unknown"

  res <- list(
    call       = object$call,
    model.info = list(
      family = "MGGUM",
      description = .fc_model_name("MGGUM", object),
      D      = D,
      method = object$method
    ),
    data.info = list(
      N = N,
      I = I,
      response.type = paste0("polytomous (", cat_range, " categories)")
    ),
    fit.stats = list(
      LogLik = logLik.val,
      npar   = object$npar,
      AIC    = AIC,
      BIC    = BIC
    ),
    par.summary   = .fc_summary_par_range(object$par$est, object$par$free, digits),
    theta.summary = .fc_summary_theta(object$theta$est, digits),
    Corr          = .fc_summary_corr(object$Corr$est, digits),
    convergence   = .fc_summary_convergence(object),
    digits        = digits
  )
  class(res) <- "summary.MGGUM"
  invisible(res)
}

# ---- summary.FCMIRT ----

#' @describeIn summary Summary method for \code{FCMIRT} objects.
#'   Forced-Choice Multidimensional IRT: dominance model with sequential
#'   ranking over item endorsement logits at the block level. Extracts
#'   block configuration, FC type
#'   (RANK/MOLE/PICK), and statement-level parameter summaries.
#' @method summary FCMIRT
#' @export
summary.FCMIRT <- function(object, digits = 4, ...) {
  N <- nrow(object$response)
  B <- length(object$block.items)
  I <- length(unlist(object$block.items))
  D <- object$arguments$D
  model <- if (!is.null(object$arguments$model)) object$arguments$model
           else if (!is.null(object$model)) object$model else "m2pl"
  logLik.val <- as.numeric(object$logLik)

  AIC <- -2 * logLik.val + 2 * object$npar
  BIC <- -2 * logLik.val + object$npar * log(N)

  fc.types <- if (!is.null(object$fc.type)) {
    paste(unique(object$fc.type), collapse = ", ")
  } else "RANK"

  res <- list(
    call       = object$call,
    model.info = list(
      family = "FCMIRT",
      description = .fc_model_name("FCMIRT", object),
      model  = model,
      D      = D,
      method = object$method,
      fc.type = fc.types
    ),
    data.info = list(
      N = N,
      B = B,
      I = I,
      response.type = paste0("forced-choice (", fc.types, ")")
    ),
    fit.stats = list(
      LogLik = logLik.val,
      npar   = object$npar,
      AIC    = AIC,
      BIC    = BIC
    ),
    par.summary   = .fc_summary_par_range(object$par$est, object$par$free, digits),
    theta.summary = .fc_summary_theta(object$theta$est, digits),
    Corr          = .fc_summary_corr(object$Corr$est, digits),
    convergence   = .fc_summary_convergence(object),
    digits        = digits
  )
  class(res) <- "summary.FCMIRT"
  invisible(res)
}

# ---- summary.FCDCM ----

#' @describeIn summary Summary method for \code{FCDCM} objects.
#'   Forced-Choice Diagnostic Classification Model: higher-order trait with
#'   exact marginalization over \eqn{2^D} attribute profiles. Extracts
#'   DCM type (DINA/DINO), higher-order delta parameters, block eta
#'   parameters, and posterior classification diagnostics.
#' @method summary FCDCM
#' @export
summary.FCDCM <- function(object, digits = 4, ...) {
  response <- object$response
  N <- nrow(response)
  B <- ncol(response)
  D <- object$arguments$D
  logLik.val <- as.numeric(object$logLik)

  AIC <- -2 * logLik.val + 2 * object$npar
  BIC <- -2 * logLik.val + object$npar * log(N)

  dcm.type <- if (!is.null(object$dcm.type)) {
    paste(unique(object$dcm.type), collapse = ", ")
  } else "DINA"

  # Posterior classification summary
  alpha.est <- object$alpha
  if (!is.null(alpha.est)) {
    if (is.list(alpha.est) && !is.null(alpha.est$est)) {
      alpha.mat <- as.matrix(alpha.est$est)
    } else if (is.list(alpha.est) && !is.null(alpha.est$prob)) {
      alpha.mat <- as.matrix(alpha.est$prob)
    } else {
      alpha.mat <- as.matrix(alpha.est)
    }
    n.patterns <- nrow(unique(alpha.mat > 0.5))
    max.post <- if (!is.null(object$class.post)) {
      mean(apply(as.matrix(object$class.post), 1, max), na.rm = TRUE)
    } else NA_real_

    class.info <- list(
      n_patterns_observed = n.patterns,
      n_patterns_possible = 2^D,
      mean_max_posterior = max.post
    )
  } else {
    class.info <- list(n_patterns_observed = NA_integer_,
                       n_patterns_possible = 2^D,
                       mean_max_posterior = NA_real_)
  }

  res <- list(
    call       = object$call,
    model.info = list(
      family = "FCDCM",
      description = .fc_model_name("FCDCM", object),
      D      = D,
      model  = dcm.type,
      method = object$method
    ),
    data.info = list(
      N = N,
      B = B,
      I = nrow(object$Q.matrix),
      response.type = "forced-choice binary blocks"
    ),
    fit.stats = list(
      LogLik = logLik.val,
      npar   = object$npar,
      AIC    = AIC,
      BIC    = BIC,
      Corr   = 1
    ),
    delta.summary = if (!is.null(object$delta$est)) {
      round(as.matrix(object$delta$est), digits)
    } else NULL,
    par.summary   = .fc_summary_par_range(object$par$est, NULL, digits),
    theta.summary = .fc_summary_theta(object$theta$est, digits),
    classification = class.info,
    convergence   = .fc_summary_convergence(object),
    digits        = digits
  )
  class(res) <- "summary.FCDCM"
  invisible(res)
}

# ---- summary.FCGDINA ----

#' @describeIn summary Summary method for \code{FCGDINA} objects.
#'   Forced-Choice GDINA model: CDM item model (DINA/DINO/ACDM/GDINA) with
#'   sequential forced-choice ranking over item endorsement logits.
#' @method summary FCGDINA
#' @export
summary.FCGDINA <- function(object, digits = 4, ...) {
  response <- object$response
  N <- nrow(response)
  B <- length(object$block.items)
  I <- length(unlist(object$block.items, use.names = FALSE))
  D <- ncol(object$Q.matrix)
  logLik.val <- as.numeric(object$logLik)

  AIC <- -2 * logLik.val + 2 * object$npar
  BIC <- -2 * logLik.val + object$npar * log(N)

  fc.types <- if (!is.null(object$fc.type)) {
    paste(unique(object$fc.type), collapse = ", ")
  } else "RANK"

  max.delta.len <- max(vapply(object$delta$est, length, integer(1L)))
  delta.mat <- do.call(rbind, lapply(object$delta$est, function(x) {
    c(x, rep(NA_real_, max.delta.len - length(x)))
  }))
  rownames(delta.mat) <- paste0("item", seq_len(nrow(delta.mat)))
  colnames(delta.mat) <- paste0("delta", seq_len(ncol(delta.mat)))

  alpha.mat <- as.matrix(object$alpha$est)
  class.info <- list(
    n_patterns_observed = nrow(unique(alpha.mat > 0.5)),
    n_patterns_possible = 2^D,
    mean_max_posterior = if (!is.null(object$class.post)) {
      mean(apply(object$class.post, 1L, max))
    } else NA_real_
  )

  pi.stats <- if (!is.null(object$pi)) {
    nz <- sum(object$pi > 1e-8)
    list(n_nonzero = nz, max = max(object$pi), min = min(object$pi),
         entropy = -sum(object$pi[object$pi > 0] * log(object$pi[object$pi > 0])))
  } else NULL

  res <- list(
    call = object$call,
    model.info = list(
      family = "FCGDINA",
      description = .fc_model_name("FCGDINA", object),
      D = D,
      model = object$model,
      method = object$method,
      fc.type = fc.types
    ),
    data.info = list(
      N = N,
      B = B,
      I = I,
      response.type = paste0("forced-choice CDM (", fc.types, ")")
    ),
    fit.stats = list(
      LogLik = logLik.val,
      npar = object$npar,
      AIC = AIC,
      BIC = BIC
    ),
    delta.summary = round(delta.mat, digits),
    alpha.summary = .fc_summary_theta(alpha.mat, digits),
    classification = class.info,
    pi.stats = pi.stats,
    convergence = .fc_summary_convergence(object),
    digits = digits
  )
  class(res) <- "summary.FCGDINA"
  invisible(res)
}

# ---- summary.FCGGUM ----

#' @describeIn summary Summary method for \code{FCGGUM} objects.
#'   Forced-Choice Generalized Graded Unfolding Model: ideal-point model
#'   with sequential ranking over binary GGUM endorsement logits. Extracts
#'   forced-choice block structure, unfolding parameters (a, delta, tau),
#'   and convergence info.
#' @method summary FCGGUM
#' @export
summary.FCGGUM <- function(object, digits = 4, ...) {
  N <- nrow(object$response)
  B <- length(object$block.items)
  I <- length(unlist(object$block.items))
  D <- object$arguments$D
  logLik.val <- as.numeric(object$logLik)

  AIC <- -2 * logLik.val + 2 * object$npar
  BIC <- -2 * logLik.val + object$npar * log(N)

  fc.types <- if (!is.null(object$fc.type)) {
    paste(unique(object$fc.type), collapse = ", ")
  } else "RANK"

  cat_range <- if (!is.null(object$length.poly)) {
    paste0(min(object$length.poly), "--", max(object$length.poly))
  } else "unknown"

  res <- list(
    call       = object$call,
    model.info = list(
      family = "FCGGUM",
      description = .fc_model_name("FCGGUM", object),
      D      = D,
      method = object$method,
      fc.type = fc.types
    ),
    data.info = list(
      N = N,
      B = B,
      I = I,
      response.type = paste0("forced-choice unfolding (", fc.types, ", ",
                             cat_range, " categories)")
    ),
    fit.stats = list(
      LogLik = logLik.val,
      npar   = object$npar,
      AIC    = AIC,
      BIC    = BIC
    ),
    par.summary   = .fc_summary_par_range(object$par$est, object$par$free, digits),
    theta.summary = .fc_summary_theta(object$theta$est, digits),
    Corr          = .fc_summary_corr(object$Corr$est, digits),
    convergence   = .fc_summary_convergence(object),
    digits        = digits
  )
  class(res) <- "summary.FCGGUM"
  invisible(res)
}

# ---- summary.TIRT ----

#' @describeIn summary Summary method for \code{TIRT} objects.
#'   Thurstonian IRT for Forced-Choice: pairwise probit comparison of latent
#'   utility differences. Extracts statement loadings (lambda), uniquenesses
#'   (psi2), pairwise gamma matrix, and factor correlations.
#' @method summary TIRT
#' @export
summary.TIRT <- function(object, digits = 4, ...) {
  N <- nrow(object$response)
  I.pairs <- ncol(object$response)
  N.block <- length(object$block.items)
  I.states <- nrow(object$par$est)
  D <- ncol(object$par$est) - 1L
  logLik.val <- as.numeric(object$logLik)

  AIC <- -2 * logLik.val + 2 * object$npar
  BIC <- -2 * logLik.val + object$npar * log(N)

  fc.types <- if (!is.null(object$fc.type)) {
    paste(unique(object$fc.type), collapse = ", ")
  } else "RANK"

  res <- list(
    call       = object$call,
    model.info = list(
      family = "TIRT",
      description = .fc_model_name("TIRT", object),
      D      = D,
      I.states = I.states,
      N.block = N.block,
      method = object$method,
      fc.type = fc.types
    ),
    data.info = list(
      N = N,
      I.pairs = I.pairs,
      I.states = I.states,
      response.type = paste0("pairwise binary (", fc.types, ")")
    ),
    fit.stats = list(
      LogLik = logLik.val,
      npar   = object$npar,
      AIC    = AIC,
      BIC    = BIC
    ),
    lambda.summary = .fc_summary_par_range(
      object$par$est[, seq_len(D), drop = FALSE], NULL, digits),
    psi2.summary = if (D + 1L <= ncol(object$par$est)) {
      matrix(summary(as.vector(object$par$est[, D + 1L])), nrow = 1L)
    } else NULL,
    theta.summary = .fc_summary_theta(object$theta$est, digits),
    Corr          = .fc_summary_corr(object$Corr$est, digits),
    convergence   = .fc_summary_convergence(object),
    digits        = digits
  )
  class(res) <- "summary.TIRT"
  invisible(res)
}


# ===========================================================================
# Internal gof table builders (used by summary.good.of.fit)
# ===========================================================================

.gof_summary_tables <- function(object, digits = 4) {
  rmsea.ci <- .gof_rmsea_ci_text(object, digits)

  tables <- list(
    Overview                  = .gof_overview_table(object, digits),
    Absolute_Fit_M2           = .gof_absolute_m2_table(object, rmsea.ci, digits),
    Comparative_Fit           = .gof_comparative_fit_table(object, digits),
    Likelihood_and_IC         = .gof_likelihood_ic_table(object, digits),
    Pseudo_R2                 = .gof_pseudo_r2_table(object, digits),
    Local_Dependence          = .gof_ld_table(object, digits),
    Classification            = .gof_class_table(object, digits)
  )

  tables[vapply(tables, function(x) nrow(x) > 0L, logical(1L))]
}

.gof_overview_table <- function(object, digits) {
  null.desc <- "independence model (item margins only)"
  if (!is.null(object$response.type) && object$response.type == "nominal") {
    null.desc <- "independence model (block margins only)"
  }
  data.frame(
    Field = c("Sample size", "Observed responses", "Items / indicators",
              "Response type", "Free parameters", "Moments",
              "Bivariate moments", "Null model"),
    Value = c(
      .gof_format(object$N, digits, int = TRUE),
      .gof_format(object$observed.responses.n, digits, int = TRUE),
      .gof_format(object$nitems, digits, int = TRUE),
      if (!is.null(object$response.type)) object$response.type else "NA",
      .gof_format(object$npar, digits, int = TRUE),
      .gof_format(object$moments.n, digits, int = TRUE),
      .gof_format(object$pairs.n, digits, int = TRUE),
      null.desc
    ),
    stringsAsFactors = FALSE
  )
}

.gof_absolute_m2_table <- function(object, rmsea.ci, digits) {
  .gof_metric_table(
    c("M2", "df", "p-value",
      "RMSEA", "RMSEA CI", "RMSEA p-close", "SRMSR",
      "McDonald NCI"),
    list(object$M2.value, object$df,
         object$p.value, object$RMSEA2, rmsea.ci,
         object$RMSEA.p.close, object$SRMSR, object$McDonald.NCI),
    c("lower", "positive", "> .05",
      "<= .05", "narrow", "> .05", "<= .08", ">= .90"),
    c("limited-information test",
      "degrees of freedom", .gof_interpret("p.value", object$p.value, object),
      .gof_interpret("RMSEA2", object$RMSEA2, object),
      "RMSEA confidence interval",
      .gof_interpret("RMSEA.p.close", object$RMSEA.p.close, object),
      .gof_interpret("SRMSR", object$SRMSR, object),
      .gof_interpret("McDonald.NCI", object$McDonald.NCI, object)),
    digits
  )
}

.gof_comparative_fit_table <- function(object, digits) {
  .gof_metric_table(
    c("CFI", "TLI", "IFI"),
    list(object$CFI, object$TLI, object$IFI),
    c(">= .90", ">= .90", ">= .90"),
    c(.gof_interpret("CFI", object$CFI, object),
      .gof_interpret("TLI", object$TLI, object),
      .gof_interpret("IFI", object$IFI, object)),
    digits
  )
}

.gof_likelihood_ic_table <- function(object, digits) {
  .gof_metric_table(
    c("SUB:Likelihood",
      "Log-likelihood", "Deviance (-2LL)",
      "SUB:Null-model Improvement",
      "G2 (vs null)", "df", "p",
      "SUB:Information Criteria",
      "AIC", "BIC", "SABIC", "CAIC"),
    list(NA, object$LogLik, object$deviance,
         NA, object$LRT, object$LRT.df, object$LRT.p,
         NA, object$AIC, object$BIC, object$SABIC, object$CAIC),
    c("", "higher", "lower",
      "", "", "", "< .05",
      "", "lower", "lower", "lower", "lower"),
    c("",
      "fitted log-likelihood", "fitted model deviance",
      "",
      "G2 = 2(LL_fitted - LL_null)", "npar (free parameters)",
      "significant => model improves over null",
      "",
      "Akaike IC", "Bayesian IC", "sample-adjusted BIC",
      "consistent AIC"),
    digits
  )
}

.gof_pseudo_r2_table <- function(object, digits) {
  .gof_metric_table(
    c("McFadden R2", "McFadden adj. R2",
      "Cox-Snell R2", "Nagelkerke R2",
      "Aldrich-Nelson R2", "Veall-Zimmermann R2"),
    list(object$McFadden.R2, object$McFadden.adj.R2,
         object$CoxSnell.R2, object$Nagelkerke.R2,
         object$AldrichNelson.R2, object$VeallZimmermann.R2),
    c("higher", "higher", "higher", "higher", "higher", "higher"),
    c("1 - (LL_fitted / LL_null)",
      "McFadden R2 penalised by npar",
      "1 - exp(-G2_null / N)",
      "(Cox-Snell) / (1 - exp(2*LL_null/N))",
      "G2_null / (G2_null + N)",
      "Aldrich-Nelson with upper-bound correction"),
    digits
  )
}

.gof_ld_table <- function(object, digits) {
  .gof_metric_table(
    c("Q3 mean", "Q3 |max|", "Adjusted Q3 |max|", "Q3 P95"),
    list(object$q3.mean, object$q3.abs.max, object$q3.adjusted.abs.max,
         object$q3.p95),
    c("near 0", "lower", "lower", "lower"),
    c("mean residual correlation", "largest absolute Q3",
      "mean-adjusted largest absolute Q3", "95th percentile Q3"),
    digits
  )
}

.gof_class_table <- function(object, digits) {
  .gof_metric_table(
    c("Normalized entropy", "Mean max posterior", "Min max posterior"),
    list(object$posterior.entropy.normalized, object$mean.max.posterior,
         object$min.max.posterior),
    c(">= .70", "higher", "higher"),
    c("classification separation", "average maximum posterior",
      "weakest case maximum posterior"),
    digits
  )
}

# ---- gof helpers ----

.gof_metric_table <- function(metric, estimate, recommended, interpretation, digits) {
  if (!is.list(estimate)) {
    estimate <- as.list(estimate)
  }
  estimate_chr <- vapply(estimate, .gof_format, character(1L), digits = digits)

  out <- data.frame(
    Metric         = metric,
    Estimate       = estimate_chr,
    Recommended    = recommended,
    Interpretation = interpretation,
    stringsAsFactors = FALSE
  )
  out
}

.gof_rmsea_ci_text <- function(object, digits) {
  lo <- object$RMSEA.lower
  hi <- object$RMSEA.upper
  if (is.null(lo) || is.null(hi) || is.na(lo) || is.na(hi)) {
    return(NA_character_)
  }
  level <- object$RMSEA.CI.level
  if (is.null(level) || is.na(level)) {
    level <- 0.90
  }
  paste0(.gof_format(100 * level, 0), "% [",
         .gof_format(lo, digits), ", ", .gof_format(hi, digits), "]")
}

.gof_format <- function(x, digits = 4, int = FALSE) {
  if (is.null(x) || length(x) == 0L) {
    return("NA")
  }
  x <- x[[1L]]
  if (is.character(x)) {
    if (is.na(x) || !nzchar(x)) {
      return("NA")
    }
    numeric.pattern <- "^[-+]?((\\d+\\.?\\d*)|(\\.\\d+))([eE][-+]?\\d+)?$"
    if (grepl(numeric.pattern, x)) {
      x <- as.numeric(x)
    } else {
      return(x)
    }
  }
  if (is.logical(x)) {
    return(as.character(x))
  }
  if (is.na(x)) {
    return("NA")
  }
  if (!is.finite(x)) {
    return(as.character(x))
  }
  if (int || abs(x - round(x)) < .Machine$double.eps^0.5) {
    return(formatC(round(x), format = "d", big.mark = ","))
  }
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

.gof_interpret <- function(metric, value, object) {
  if (length(value) == 0L || is.null(value) || is.na(value)) {
    return("not available")
  }
  switch(metric,
    p.value = if (value >= 0.05) "limited-information misfit not significant"
              else "significant limited-information misfit",
    M2.value = "limited-information marginal discrepancy",
    df = "degrees of freedom for the fit statistic",
    RMSEA2 = if (value <= 0.05) "close approximate fit"
             else if (value <= 0.08) "reasonable approximate fit"
             else if (value <= 0.10) "mediocre approximate fit"
             else "poor approximate fit",
    RMSEA.CI = "interval estimate for approximate-fit error",
    RMSEA.p.close = if (value >= 0.05) "close-fit hypothesis not rejected"
                    else "close-fit hypothesis rejected",
    SRMSR = if (value <= 0.08) "small standardized marginal residuals"
            else "large standardized marginal residuals",
    McDonald.NCI = if (value >= 0.95) "very high noncentrality fit"
                   else if (value >= 0.90) "acceptable noncentrality fit"
                   else "low noncentrality fit",
    CFI = .gof_interpret_incremental(value),
    TLI = .gof_interpret_incremental(value),
    NFI = .gof_interpret_incremental(value),
    IFI = .gof_interpret_incremental(value),
    RFI = .gof_interpret_incremental(value),
    pattern.p = if (value >= 0.05) "full-pattern misfit not significant"
                else "significant full-pattern misfit",
    "not available"
  )
}

.gof_interpret_incremental <- function(value) {
  if (is.na(value)) {
    return("not available")
  }
  if (value >= 0.95) {
    "good relative fit by common heuristic"
  } else if (value >= 0.90) {
    "acceptable relative fit by common heuristic"
  } else {
    "weak relative fit by common heuristic"
  }
}
