#' @title S3 Methods: vcov
#'
#' @description
#' Returns the variance--covariance matrix of the estimated model
#' parameters from a fitted \pkg{ForceChoice} model object. The package stores
#' marginal standard errors rather than a joint covariance estimate, so the
#' returned matrix is \eqn{\mathrm{diag}(\mathrm{SE}^2)} for every backend.
#'
#' @param object A fitted model object of class \code{"MIRT"},
#'   \code{"MGPCM"}, \code{"MGGUM"}, \code{"FCMIRT"}, \code{"FCDCM"},
#'   \code{"FCGDINA"}, \code{"FCGGUM"}, or \code{"TIRT"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A square numeric matrix. Rows and columns correspond to the
#'   free model parameters in the package's canonical parameter map.
#'
#' @details
#' Standard errors reflect the active estimation backend. Off-diagonal
#' covariances are not retained in fitted objects and are therefore reported
#' as zero. Parameters whose standard errors were not estimated have
#' \code{NA} diagonal entries.
#'
#' @name vcov
NULL


#' @title S3 Methods: confint
#'
#' @description
#' Computes Wald-type confidence intervals for the model parameters of a
#' fitted \pkg{ForceChoice} model object. Intervals are constructed as
#' \eqn{\hat{\psi}_k \pm z_{1-\alpha/2}\,\mathrm{SE}_k}, where
#' \eqn{\hat{\psi}_k} is the point estimate and \eqn{\mathrm{SE}_k} is
#' the standard error of the \eqn{k}-th free parameter.
#'
#' @param object A fitted model object.
#' @param parm Optional numeric indices or character names selecting free
#'   parameters. By default, intervals are returned for all free parameters.
#' @param level Numeric; confidence level (default: \code{0.95} for 95%
#'   intervals).
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric matrix with two columns, \code{2.5%} and
#'   \code{97.5%} (column names adjust automatically for non-default
#'   \code{level}). Each row corresponds to one free parameter.
#'
#' @name confint
NULL


#' @title S3 Methods: nobs
#'
#' @description
#' Extracts the number of observations (persons) from a fitted
#' \pkg{ForceChoice} model object. Returns the number of rows in the
#' response matrix used for model fitting.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An integer: the sample size \eqn{N}.
#'
#' @name nobs
#' @importFrom stats nobs
NULL


#' @title S3 Methods: deviance
#'
#' @description
#' Returns the model deviance, defined as \eqn{-2 \ell} where \eqn{\ell}
#' is the marginal log-likelihood evaluated at the parameter estimates.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric scalar: \eqn{-2 \log L}.
#'
#' @name deviance
NULL


# ---- Internal helpers ----

.fc_matrix_parameter_map <- function(est, se, free = NULL, prefix = "par") {
  est <- as.matrix(est)
  if (is.null(se)) {
    se <- matrix(NA_real_, nrow(est), ncol(est), dimnames = dimnames(est))
  } else {
    se <- as.matrix(se)
  }
  if (!identical(dim(est), dim(se))) {
    stop("Estimate and standard-error matrices have different dimensions.",
         call. = FALSE)
  }
  if (is.null(free)) free <- matrix(TRUE, nrow(est), ncol(est))
  free <- as.matrix(free)
  if (!identical(dim(est), dim(free))) {
    stop("Estimate and free-parameter matrices have different dimensions.",
         call. = FALSE)
  }

  idx <- which(free, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    return(data.frame(name = character(), estimate = numeric(), se = numeric()))
  }
  row_labels <- rownames(est)
  if (is.null(row_labels)) row_labels <- as.character(seq_len(nrow(est)))
  col_labels <- colnames(est)
  if (is.null(col_labels)) col_labels <- as.character(seq_len(ncol(est)))
  data.frame(
    name = paste0(prefix, ".", col_labels[idx[, 2L]], "[",
                  row_labels[idx[, 1L]], "]"),
    estimate = as.numeric(est[free]),
    se = as.numeric(se[free]),
    stringsAsFactors = FALSE
  )
}

.fc_vector_parameter_map <- function(est, se = NULL, free = NULL,
                                     prefix) {
  labels <- names(est)
  est <- as.numeric(est)
  if (is.null(se)) se <- rep(NA_real_, length(est))
  se <- as.numeric(se)
  if (length(se) != length(est)) {
    stop("Estimate and standard-error vectors have different lengths.",
         call. = FALSE)
  }
  if (is.null(free)) free <- rep(TRUE, length(est))
  if (length(free) != length(est)) {
    stop("Estimate and free-parameter vectors have different lengths.",
         call. = FALSE)
  }
  if (is.null(labels)) labels <- as.character(seq_along(est))
  data.frame(
    name = paste0(prefix, "[", labels[free], "]"),
    estimate = est[free],
    se = se[free],
    stringsAsFactors = FALSE
  )
}

.fc_corr_parameter_map <- function(object) {
  Corr <- object$Corr$est
  if (is.null(Corr) || nrow(as.matrix(Corr)) < 2L) {
    return(data.frame(name = character(), estimate = numeric(), se = numeric()))
  }
  Corr <- as.matrix(Corr)
  Corr.se <- object$Corr$se
  if (is.null(Corr.se)) {
    Corr.se <- matrix(NA_real_, nrow(Corr), ncol(Corr),
                      dimnames = dimnames(Corr))
  }
  .fc_matrix_parameter_map(
    Corr, Corr.se, free = lower.tri(Corr), prefix = "Corr"
  )
}

.fc_fcgdina_stan_delta_free_map <- function(object) {
  rank <- object$delta.identification$rank
  if (!is.numeric(rank) || length(rank) != 1L ||
      !is.finite(rank) || rank < 1L) {
    stop("FCGDINA Stan object does not contain a valid delta rank.",
         call. = FALSE)
  }
  rank <- as.integer(rank)

  draws <- object$MCMC.obj$delta_free
  if (!is.null(draws)) {
    draws <- as.matrix(draws)
    if (ncol(draws) != rank && nrow(draws) == rank) {
      draws <- t(draws)
    }
    if (ncol(draws) != rank) {
      stop("FCGDINA Stan 'delta_free' draws are inconsistent with rank.",
           call. = FALSE)
    }
    est <- colMeans(draws)
    se <- if (nrow(draws) > 1L) apply(draws, 2L, stats::sd) else
      rep(NA_real_, rank)
  } else if (!is.null(object$stan.obj)) {
    stan.sum <- rstan::summary(object$stan.obj, pars = "delta_free")$summary
    if (nrow(stan.sum) != rank) {
      stop("FCGDINA Stan 'delta_free' summary is inconsistent with rank.",
           call. = FALSE)
    }
    est <- stan.sum[, "mean"]
    se <- stan.sum[, "sd"]
  } else {
    stop("FCGDINA Stan object does not contain delta_free draws.",
         call. = FALSE)
  }

  data.frame(
    name = paste0("delta_free[", seq_len(rank), "]"),
    estimate = as.numeric(est),
    se = as.numeric(se),
    stringsAsFactors = FALSE
  )
}

.fc_parameter_map <- function(object) {
  cls <- class(object)[1L]
  if (is.null(object$npar) || length(object$npar) != 1L ||
      !is.finite(object$npar)) {
    stop("Fitted object does not contain a valid 'npar'.", call. = FALSE)
  }

  map <- switch(cls,
    MIRT =, MGPCM =, MGGUM =, FCMIRT =, FCGGUM = {
      .fc_matrix_parameter_map(
        object$par$est, object$par$se, object$par$free, prefix = "par"
      )
    },
    FCDCM = {
      rbind(
        .fc_matrix_parameter_map(
          object$delta$est, object$delta$se, prefix = "delta"
        ),
        .fc_matrix_parameter_map(
          object$par$est, object$par$se, object$par$free, prefix = "par"
        )
      )
    },
    FCGDINA = {
      delta_map <- if (!is.null(object$delta.identification)) {
        .fc_fcgdina_stan_delta_free_map(object)
      } else {
        delta_maps <- Map(
          function(est, se, i) {
            .fc_vector_parameter_map(est, se, prefix = paste0("delta", i))
          },
          object$delta$est, object$delta$se, seq_along(object$delta$est)
        )
        do.call(rbind, delta_maps)
      }
      C <- length(object$pi)
      pi_map <- .fc_vector_parameter_map(
        object$pi, object$pi.se, free = seq_len(C) < C, prefix = "pi"
      )
      rbind(delta_map, pi_map)
    },
    TIRT = {
      rbind(
        .fc_matrix_parameter_map(
          object$par$est, object$par$se, object$par$free, prefix = "par"
        ),
        .fc_vector_parameter_map(
          object$gamma$est, object$gamma$se, object$gamma$free,
          prefix = "gamma"
        )
      )
    },
    stop("Unsupported fitted-object class: ", cls, call. = FALSE)
  )

  remaining <- as.integer(object$npar) - nrow(map)
  if (remaining > 0L && cls %in% c(
    "MIRT", "MGPCM", "MGGUM", "FCMIRT", "FCGGUM", "TIRT"
  )) {
    corr_map <- .fc_corr_parameter_map(object)
    if (nrow(corr_map) != remaining) {
      stop("Free correlation-parameter count is inconsistent with 'npar'.",
           call. = FALSE)
    }
    map <- rbind(map, corr_map)
  }
  if (nrow(map) != as.integer(object$npar)) {
    stop("Canonical free-parameter map is inconsistent with 'npar'.",
         call. = FALSE)
  }
  rownames(map) <- NULL
  map
}

.fc_vcov_impl <- function(object, ...) {
  map <- .fc_parameter_map(object)
  out <- diag(map$se^2, nrow = nrow(map))
  dimnames(out) <- list(map$name, map$name)
  out
}

.fc_confint_impl <- function(object, parm = NULL, level = 0.95, ...) {
  if (length(level) != 1L || !is.numeric(level) || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a numeric scalar strictly between 0 and 1.",
         call. = FALSE)
  }
  map <- .fc_parameter_map(object)
  if (!is.null(parm)) {
    if (is.character(parm)) {
      idx <- match(parm, map$name)
      if (anyNA(idx)) {
        stop("Unknown parameter name(s): ",
             paste(parm[is.na(idx)], collapse = ", "), call. = FALSE)
      }
    } else if (is.numeric(parm) && all(is.finite(parm)) &&
               all(parm == floor(parm)) && all(parm >= 1L) &&
               all(parm <= nrow(map))) {
      idx <- as.integer(parm)
    } else {
      stop("'parm' must contain valid parameter names or indices.",
           call. = FALSE)
    }
    map <- map[idx, , drop = FALSE]
  }
  alpha <- (1 - level) / 2
  z <- stats::qnorm(1 - alpha)
  probs <- c(alpha, 1 - alpha)
  out <- cbind(map$estimate - z * map$se, map$estimate + z * map$se)
  colnames(out) <- paste0(format(100 * probs, trim = TRUE), " %")
  rownames(out) <- map$name
  out
}

.fc_nobs_impl <- function(object, ...) {
  response <- object$response
  if (is.null(response)) response <- object$arguments$data
  if (is.null(response)) response <- object$arguments$response
  if (is.null(response)) {
    stop("Fitted object does not contain response data.", call. = FALSE)
  }
  nrow(response)
}

.fc_deviance_impl <- function(object, ...) {
  value <- object$logLik
  if (is.null(value)) value <- logLik(object)
  -2 * as.numeric(value)
}


# ===========================================================================
# vcov
# ===========================================================================

#' @describeIn vcov Variance--covariance matrix of parameter estimates
#'   for \code{MIRT} objects.
#' @method vcov MIRT
#' @export
vcov.MIRT    <- .fc_vcov_impl

#' @describeIn vcov Variance--covariance matrix for \code{MGPCM} objects.
#' @method vcov MGPCM
#' @export
vcov.MGPCM   <- .fc_vcov_impl

#' @describeIn vcov Variance--covariance matrix for \code{MGGUM} objects.
#' @method vcov MGGUM
#' @export
vcov.MGGUM   <- .fc_vcov_impl

#' @describeIn vcov Variance--covariance matrix for \code{FCMIRT} objects.
#' @method vcov FCMIRT
#' @export
vcov.FCMIRT  <- .fc_vcov_impl

#' @describeIn vcov Variance--covariance matrix for \code{FCDCM} objects.
#' @method vcov FCDCM
#' @export
vcov.FCDCM   <- .fc_vcov_impl

#' @describeIn vcov Variance--covariance matrix for \code{FCGDINA} objects.
#' @method vcov FCGDINA
#' @export
vcov.FCGDINA <- .fc_vcov_impl

#' @describeIn vcov Variance--covariance matrix for \code{FCGGUM} objects.
#' @method vcov FCGGUM
#' @export
vcov.FCGGUM  <- .fc_vcov_impl

#' @describeIn vcov Variance--covariance matrix for \code{TIRT} objects.
#' @method vcov TIRT
#' @export
vcov.TIRT    <- .fc_vcov_impl


# ===========================================================================
# confint
# ===========================================================================

#' @describeIn confint Confidence intervals for \code{MIRT} parameters.
#' @method confint MIRT
#' @export
confint.MIRT   <- .fc_confint_impl

#' @describeIn confint Confidence intervals for \code{MGPCM} parameters.
#' @method confint MGPCM
#' @export
confint.MGPCM  <- .fc_confint_impl

#' @describeIn confint Confidence intervals for \code{MGGUM} parameters.
#' @method confint MGGUM
#' @export
confint.MGGUM  <- .fc_confint_impl

#' @describeIn confint Confidence intervals for \code{FCMIRT} parameters.
#' @method confint FCMIRT
#' @export
confint.FCMIRT <- .fc_confint_impl

#' @describeIn confint Confidence intervals for \code{FCDCM} parameters.
#' @method confint FCDCM
#' @export
confint.FCDCM  <- .fc_confint_impl

#' @describeIn confint Confidence intervals for \code{FCGDINA} parameters.
#' @method confint FCGDINA
#' @export
confint.FCGDINA <- .fc_confint_impl

#' @describeIn confint Confidence intervals for \code{FCGGUM} parameters.
#' @method confint FCGGUM
#' @export
confint.FCGGUM <- .fc_confint_impl

#' @describeIn confint Confidence intervals for \code{TIRT} parameters.
#' @method confint TIRT
#' @export
confint.TIRT   <- .fc_confint_impl


# ===========================================================================
# nobs
# ===========================================================================

#' @describeIn nobs Number of observations (persons) for \code{MIRT} objects.
#' @method nobs MIRT
#' @export
nobs.MIRT    <- .fc_nobs_impl

#' @describeIn nobs Number of observations for \code{MGPCM} objects.
#' @method nobs MGPCM
#' @export
nobs.MGPCM   <- .fc_nobs_impl

#' @describeIn nobs Number of observations for \code{MGGUM} objects.
#' @method nobs MGGUM
#' @export
nobs.MGGUM   <- .fc_nobs_impl

#' @describeIn nobs Number of observations for \code{FCMIRT} objects.
#' @method nobs FCMIRT
#' @export
nobs.FCMIRT  <- .fc_nobs_impl

#' @describeIn nobs Number of observations for \code{FCDCM} objects.
#' @method nobs FCDCM
#' @export
nobs.FCDCM   <- .fc_nobs_impl

#' @describeIn nobs Number of observations for \code{FCGDINA} objects.
#' @method nobs FCGDINA
#' @export
nobs.FCGDINA <- .fc_nobs_impl

#' @describeIn nobs Number of observations for \code{FCGGUM} objects.
#' @method nobs FCGGUM
#' @export
nobs.FCGGUM  <- .fc_nobs_impl

#' @describeIn nobs Number of observations for \code{TIRT} objects.
#' @method nobs TIRT
#' @export
nobs.TIRT    <- .fc_nobs_impl


# ===========================================================================
# deviance
# ===========================================================================

#' @describeIn deviance Model deviance (\eqn{-2 \ell}) for \code{MIRT} objects.
#' @method deviance MIRT
#' @export
deviance.MIRT   <- .fc_deviance_impl

#' @describeIn deviance Model deviance for \code{MGPCM} objects.
#' @method deviance MGPCM
#' @export
deviance.MGPCM  <- .fc_deviance_impl

#' @describeIn deviance Model deviance for \code{MGGUM} objects.
#' @method deviance MGGUM
#' @export
deviance.MGGUM  <- .fc_deviance_impl

#' @describeIn deviance Model deviance for \code{FCMIRT} objects.
#' @method deviance FCMIRT
#' @export
deviance.FCMIRT <- .fc_deviance_impl

#' @describeIn deviance Model deviance for \code{FCDCM} objects.
#' @method deviance FCDCM
#' @export
deviance.FCDCM  <- .fc_deviance_impl

#' @describeIn deviance Model deviance for \code{FCGDINA} objects.
#' @method deviance FCGDINA
#' @export
deviance.FCGDINA <- .fc_deviance_impl

#' @describeIn deviance Model deviance for \code{FCGGUM} objects.
#' @method deviance FCGGUM
#' @export
deviance.FCGGUM <- .fc_deviance_impl

#' @describeIn deviance Model deviance for \code{TIRT} objects.
#' @method deviance TIRT
#' @export
deviance.TIRT   <- .fc_deviance_impl
