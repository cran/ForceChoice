#' @title S3 Methods: fitted
#'
#' @description
#' Extract model-implied response probabilities from fitted
#' \pkg{ForceChoice} model objects. Returns the conditional probability of
#' each observed response given each person's estimated latent trait
#' \eqn{\hat{\boldsymbol{\theta}}_j}, evaluated via the model's response
#' function.
#'
#' @param object A fitted model object of class \code{"MIRT"},
#'   \code{"MGPCM"}, \code{"MGGUM"}, \code{"FCMIRT"}, \code{"FCDCM"},
#'   \code{"FCGDINA"}, \code{"FCGGUM"}, or \code{"TIRT"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric matrix of model-implied probabilities. Dimensions:
#'   \itemize{
#'     \item Binary models (MIRT, TIRT): \eqn{N \times I}
#'     \item Polytomous models (MGPCM, MGGUM): \eqn{N \times \sum K_i}
#'           (stacked category probabilities)
#'     \item Forced-choice models (FCMIRT, FCGDINA, FCGGUM): \eqn{N \times \sum K_b}
#'           (block-pattern probabilities)
#'     \item FCDCM: \eqn{N \times B} (marginal block-choice probabilities
#'           after integrating over \eqn{2^D} attribute patterns)
#'   }
#'
#' @details
#' For binary-response models, each entry is the endorsement probability
#' \eqn{P(Y_{ij}=1 \mid \hat{\boldsymbol{\theta}}_j)}. For polytomous
#' models, probabilities are returned in stacked format (one column per
#' category). For forced-choice models, each entry is the probability of
#' the observed ranking pattern in the given block.
#'
#' @name fitted
NULL


#' @title S3 Methods: residuals
#'
#' @description
#' Compute residuals from fitted \pkg{ForceChoice} model objects.
#' Supported types are raw (response), Pearson, and deviance residuals.
#'
#' @param object A fitted model object.
#' @param type Character; \code{"raw"} (default), \code{"pearson"}, or
#'   \code{"deviance"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric matrix with the same dimensions as the corresponding
#'   \code{fitted} output.
#'
#' @details
#' Residual types follow standard GLM conventions for binary responses:
#' \itemize{
#'   \item \strong{raw}: \eqn{e_{ij} = y_{ij} - \hat{p}_{ij}}
#'   \item \strong{Pearson}:
#'         \eqn{r_{ij} = e_{ij} / \sqrt{\hat{p}_{ij}(1-\hat{p}_{ij})}}
#'   \item \strong{deviance}:
#'         \eqn{d_{ij} = \mathrm{sign}(e_{ij}) \sqrt{2\left[y_{ij}\log\frac{y_{ij}}{\hat{p}_{ij}} + (1-y_{ij})\log\frac{1-y_{ij}}{1-\hat{p}_{ij}}\right]}}
#' }
#' For polytomous and forced-choice models, residuals are computed on a
#' category-indicator matrix. This treats ordered categories and nominal
#' forced-choice response patterns as multinomial outcomes rather than as
#' numeric scores.
#'
#' @name residuals
NULL


# ---- Internal helpers ----

.fc_fitted_binary <- function(object) {
  theta <- as.matrix(object$theta$est)
  par   <- as.matrix(object$par$est)
  # Dispatch to the appropriate model.* function based on class
  cls <- class(object)[1L]
  switch(cls,
    MIRT   = model.MIRT(theta, par),
    TIRT   = {
      D <- ncol(par) - 1L
      lambda <- unlist(lapply(seq_len(nrow(par)), function(x) {
        d <- which(object$Q.matrix[x, ] != 0)[1]
        return(par[x, d])
      }))
      psi2 <- par[, D + 1L]
      model.TIRT(theta, lambda, psi2, object$gamma.matrix$est,
                 object$Q.matrix, object$pairs.matrix)
    },
    stop("fitted() not implemented for class ", cls, call. = FALSE)
  )
}

.fc_fitted_polytomous <- function(object) {
  theta <- as.matrix(object$theta$est)
  par   <- as.matrix(object$par$est)
  cls <- class(object)[1L]
  switch(cls,
    MGPCM = model.MGPCM(theta, par),
    MGGUM = model.MGGUM(theta, par),
    stop("fitted() not implemented for class ", cls, call. = FALSE)
  )
}

.fc_fitted_fc <- function(object) {
  theta <- as.matrix(object$theta$est)
  par   <- as.matrix(object$par$est)
  cls <- class(object)[1L]
  switch(cls,
    FCMIRT = model.FCMIRT(theta, par, object$patterns.total, object$patterns),
    FCGGUM = model.FCGGUM(theta, par, object$patterns.total, object$patterns),
    stop("fitted() not implemented for class ", cls, call. = FALSE)
  )
}

.fc_fitted_fcdcm <- function(object) {
  theta <- as.matrix(object$theta$est)
  # P(alpha_k | theta_j) for all j = 1..N, k = 1..2^D
  alpha_prob <- fcdcm_alpha_profile_prob(
    theta = theta,
    delta1 = object$delta$est[, "delta1"],
    delta0 = object$delta$est[, "delta0"],
    alpha.patterns = object$alpha.patterns
  )
  # Block-choice probabilities for each zeta (condensation) pattern
  prob_per_pattern <- model.FCDCM(
    zeta = object$zeta.patterns,
    par = object$par$est,
    patterns = object$patterns
  )
  # Marginalize over alpha patterns:
  #   P(Y_jb = 1 | theta_j) = sum_k P(Y_b = 1 | zeta_k) * P(alpha_k | theta_j)
  alpha_prob %*% prob_per_pattern
}

.fc_fitted_fcgdina <- function(object) {
  prep <- fcgdina_prepare(
    data = object$response,
    Q.matrix = object$Q.matrix,
    model = object$model,
    block.items = object$block.items,
    fc.type = object$fc.type
  )
  prob.class <- fcgdina_prob_class_link(
    prep, object$delta$est, link = fcgdina_delta_link(object))
  class.prob <- object$class.post
  class.prob %*% prob.class
}

.fc_residuals_raw <- function(fitted, response) {
  response - fitted
}

.fc_residuals_pearson <- function(fitted, response) {
  raw <- response - fitted
  sd  <- sqrt(pmax(fitted * (1 - fitted), .Machine$double.eps))
  raw / sd
}

.fc_residuals_deviance <- function(fitted, response) {
  eps <- .Machine$double.eps
  p <- pmax(pmin(fitted, 1 - eps), eps)
  y <- response
  dev_contrib <- 2 * (y * log(pmax(y, eps) / p) +
                      (1 - y) * log(pmax(1 - y, eps) / (1 - p)))
  sign(y - p) * sqrt(pmax(dev_contrib, 0))
}

.fc_check_residual_dims <- function(fitted, response) {
  if (!identical(dim(fitted), dim(response))) {
    stop(
      "'fitted' and expanded response matrices must have identical dimensions.",
      call. = FALSE
    )
  }
}

.fc_binary_response <- function(response, fitted) {
  response <- as.matrix(response)
  if (anyNA(response) || any(response != floor(response)) ||
      !all(response %in% c(0, 1))) {
    stop("Binary residuals require a response matrix coded 0/1.",
         call. = FALSE)
  }
  storage.mode(response) <- "numeric"
  .fc_check_residual_dims(fitted, response)
  dimnames(response) <- dimnames(fitted)
  response
}

.fc_response_indicators <- function(response, n_cat, first = 0L,
                                    prefix = "Y") {
  response <- as.matrix(response)
  if (anyNA(response) || any(response != floor(response))) {
    stop("Expanded residuals require finite integer response categories.",
         call. = FALSE)
  }
  storage.mode(response) <- "integer"

  N <- nrow(response)
  I <- ncol(response)
  n_cat <- as.integer(n_cat)
  if (length(n_cat) != I || anyNA(n_cat) || any(n_cat < 2L)) {
    stop("'n_cat' must contain one category count >= 2 per response column.",
         call. = FALSE)
  }

  first <- as.integer(first)
  if (length(first) != 1L || is.na(first)) {
    stop("'first' must be a scalar integer category origin.", call. = FALSE)
  }

  out <- matrix(0, N, sum(n_cat))
  rownames(out) <- rownames(response)
  col.names <- character(sum(n_cat))
  base.names <- colnames(response)
  if (is.null(base.names)) base.names <- paste0(prefix, seq_len(I))

  idx <- 0L
  for (i in seq_len(I)) {
    yi <- response[, i]
    last <- first + n_cat[i] - 1L
    if (any(yi < first | yi > last)) {
      stop(
        "Response column ", i, " contains categories outside ",
        first, ":", last, ".",
        call. = FALSE
      )
    }
    out[cbind(seq_len(N), idx + yi - first + 1L)] <- 1
    col.names[idx + seq_len(n_cat[i])] <-
      paste0(base.names[i], ".", first:last)
    idx <- idx + n_cat[i]
  }

  colnames(out) <- col.names
  out
}

.fc_residuals <- function(fitted, response, type) {
  .fc_check_residual_dims(fitted, response)
  switch(type,
    raw      = .fc_residuals_raw(fitted, response),
    pearson  = .fc_residuals_pearson(fitted, response),
    deviance = .fc_residuals_deviance(fitted, response)
  )
}


# ===========================================================================
# MIRT
# ===========================================================================

#' @describeIn fitted Fitted probabilities for \code{MIRT} objects.
#'   Returns \eqn{N \times I} matrix of endorsement probabilities.
#' @method fitted MIRT
#' @export
fitted.MIRT <- function(object, ...) {
  .fc_fitted_binary(object)
}

#' @describeIn residuals Residuals for \code{MIRT} objects.
#' @method residuals MIRT
#' @export
residuals.MIRT <- function(object, type = c("raw", "pearson", "deviance"), ...) {
  type <- match.arg(type)
  fit <- fitted(object)
  resp <- .fc_binary_response(object$arguments$response, fit)
  .fc_residuals(fit, resp, type)
}

# ===========================================================================
# MGPCM
# ===========================================================================

#' @describeIn fitted Fitted probabilities for \code{MGPCM} objects.
#'   Returns stacked category probability matrix.
#' @method fitted MGPCM
#' @export
fitted.MGPCM <- function(object, ...) {
  .fc_fitted_polytomous(object)
}

#' @describeIn residuals Residuals for \code{MGPCM} objects.
#' @method residuals MGPCM
#' @export
residuals.MGPCM <- function(object, type = c("raw", "pearson", "deviance"), ...) {
  type <- match.arg(type)
  fit <- fitted(object)
  resp <- .fc_response_indicators(
    object$arguments$response, object$length.poly, first = 0L, prefix = "item"
  )
  dimnames(resp) <- dimnames(fit)
  .fc_residuals(fit, resp, type)
}

# ===========================================================================
# MGGUM
# ===========================================================================

#' @describeIn fitted Fitted probabilities for \code{MGGUM} objects.
#' @method fitted MGGUM
#' @export
fitted.MGGUM <- function(object, ...) {
  .fc_fitted_polytomous(object)
}

#' @describeIn residuals Residuals for \code{MGGUM} objects.
#' @method residuals MGGUM
#' @export
residuals.MGGUM <- function(object, type = c("raw", "pearson", "deviance"), ...) {
  type <- match.arg(type)
  fit <- fitted(object)
  resp <- .fc_response_indicators(
    object$arguments$response, object$length.poly, first = 0L, prefix = "item"
  )
  dimnames(resp) <- dimnames(fit)
  .fc_residuals(fit, resp, type)
}

# ===========================================================================
# FCMIRT
# ===========================================================================

#' @describeIn fitted Fitted probabilities for \code{FCMIRT} objects.
#'   Returns block-pattern probability matrix.
#' @method fitted FCMIRT
#' @export
fitted.FCMIRT <- function(object, ...) {
  .fc_fitted_fc(object)
}

#' @describeIn residuals Residuals for \code{FCMIRT} objects.
#' @method residuals FCMIRT
#' @export
residuals.FCMIRT <- function(object, type = c("raw", "pearson", "deviance"), ...) {
  type <- match.arg(type)
  fit <- fitted(object)
  resp <- .fc_response_indicators(
    object$response,
    vapply(object$patterns, nrow, integer(1L)),
    first = 1L,
    prefix = "block"
  )
  dimnames(resp) <- dimnames(fit)
  .fc_residuals(fit, resp, type)
}

# ===========================================================================
# FCDCM
# ===========================================================================

#' @describeIn fitted Fitted probabilities for \code{FCDCM} objects.
#'   Returns \eqn{N \times B} matrix of marginal block-choice probabilities
#'   after integrating over \eqn{2^D} attribute patterns.
#' @method fitted FCDCM
#' @export
fitted.FCDCM <- function(object, ...) {
  .fc_fitted_fcdcm(object)
}

#' @describeIn residuals Residuals for \code{FCDCM} objects.
#' @method residuals FCDCM
#' @export
residuals.FCDCM <- function(object, type = c("raw", "pearson", "deviance"), ...) {
  type <- match.arg(type)
  fit <- fitted(object)
  resp <- .fc_binary_response(object$response, fit)
  .fc_residuals(fit, resp, type)
}

# ===========================================================================
# FCGDINA
# ===========================================================================

#' @describeIn fitted Fitted probabilities for \code{FCGDINA} objects.
#'   Returns block-pattern probability matrix after integrating over posterior
#'   attribute classes.
#' @method fitted FCGDINA
#' @export
fitted.FCGDINA <- function(object, ...) {
  .fc_fitted_fcgdina(object)
}

#' @describeIn residuals Residuals for \code{FCGDINA} objects.
#' @method residuals FCGDINA
#' @export
residuals.FCGDINA <- function(object, type = c("raw", "pearson", "deviance"), ...) {
  type <- match.arg(type)
  fit <- fitted(object)
  resp <- .fc_response_indicators(
    object$response,
    vapply(object$patterns, nrow, integer(1L)),
    first = 1L,
    prefix = "block"
  )
  dimnames(resp) <- dimnames(fit)
  .fc_residuals(fit, resp, type)
}

# ===========================================================================
# FCGGUM
# ===========================================================================

#' @describeIn fitted Fitted probabilities for \code{FCGGUM} objects.
#' @method fitted FCGGUM
#' @export
fitted.FCGGUM <- function(object, ...) {
  .fc_fitted_fc(object)
}

#' @describeIn residuals Residuals for \code{FCGGUM} objects.
#' @method residuals FCGGUM
#' @export
residuals.FCGGUM <- function(object, type = c("raw", "pearson", "deviance"), ...) {
  type <- match.arg(type)
  fit <- fitted(object)
  resp <- .fc_response_indicators(
    object$response,
    vapply(object$patterns, nrow, integer(1L)),
    first = 1L,
    prefix = "block"
  )
  dimnames(resp) <- dimnames(fit)
  .fc_residuals(fit, resp, type)
}

# ===========================================================================
# TIRT
# ===========================================================================

#' @describeIn fitted Fitted probabilities for \code{TIRT} objects.
#'   Returns \eqn{N \times I_{pairs}} matrix of pairwise comparison
#'   probabilities.
#' @method fitted TIRT
#' @export
fitted.TIRT <- function(object, ...) {
  .fc_fitted_binary(object)
}

#' @describeIn residuals Residuals for \code{TIRT} objects.
#' @method residuals TIRT
#' @export
residuals.TIRT <- function(object, type = c("raw", "pearson", "deviance"), ...) {
  type <- match.arg(type)
  fit <- fitted(object)
  resp <- .fc_binary_response(object$response, fit)
  .fc_residuals(fit, resp, type)
}
