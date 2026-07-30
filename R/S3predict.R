#' @title S3 Methods: predict
#'
#' @description
#' Predict response probabilities from fitted \pkg{ForceChoice} model
#' objects. Given a fitted model and optionally new person parameter
#' (latent trait) values, returns the model-implied response probabilities.
#' When \code{newdata} is not provided, predictions are made at the
#' estimated latent trait values for each person.
#'
#' @param object A fitted model object of class \code{"MIRT"},
#'   \code{"MGPCM"}, \code{"MGGUM"}, \code{"FCMIRT"}, \code{"FCDCM"},
#'   \code{"FCGDINA"}, \code{"FCGGUM"}, or \code{"TIRT"}.
#' @param newdata An optional \eqn{M \times D} matrix of latent trait
#'   values at which to compute predictions. If \code{NULL} (default),
#'   predictions use the estimated \code{object$theta$est}.
#' @param type Character; \code{"probability"} (default) returns
#'   response probabilities, \code{"response"} returns expected category
#'   or simulated binary responses (random draws).
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' \describe{
#'   \item{\code{type = "probability"}}{For binary models (MIRT, TIRT):
#'     an \eqn{M \times I} matrix of endorsement probabilities. For
#'     polytomous models (MGPCM, MGGUM): an \eqn{M \times \sum K_i}
#'     matrix of stacked category probabilities. For forced-choice
#'     models (FCMIRT, FCGGUM): an \eqn{M \times \sum K_b} matrix of
#'     block-pattern probabilities. For FCDCM: an \eqn{M \times B}
#'     matrix of marginal block-choice probabilities.}
#'   \item{\code{type = "response"}}{For binary models: an
#'     \eqn{M \times I} matrix of simulated 0/1 responses. For
#'     polytomous models: an \eqn{M \times I} matrix of category
#'     indices. For forced-choice models: an \eqn{M \times B} matrix
#'     of chosen pattern indices.}
#' }
#'
#' @name predict
NULL


# ---- Internal helpers ----

.fc_predict_prob <- function(object, newdata) {
  cls <- class(object)[1L]
  if (cls == "FCGDINA") {
    prep <- fcgdina_prepare(
      data = object$response,
      Q.matrix = object$Q.matrix,
      model = object$model,
      block.items = object$block.items,
      fc.type = object$fc.type
    )
    prob.class <- fcgdina_prob_class_link(
      prep, object$delta$est, link = fcgdina_delta_link(object))
    if (is.null(newdata)) {
      class.prob <- object$class.post
    } else {
      nd <- as.matrix(newdata)
      if (ncol(nd) == prep$D && all(nd %in% c(0, 1))) {
        key <- pattern_key(nd)
        cls.idx <- match(key, pattern_key(prep$alpha.patterns))
        if (anyNA(cls.idx)) stop("'newdata' contains invalid alpha patterns.",
                                 call. = FALSE)
        class.prob <- matrix(0, nrow(nd), prep$C)
        class.prob[cbind(seq_len(nrow(nd)), cls.idx)] <- 1
      } else if (ncol(nd) == prep$C) {
        rs <- rowSums(nd)
        if (any(!is.finite(rs)) || any(rs <= 0)) {
          stop("'newdata' class-probability rows must have positive sums.",
               call. = FALSE)
        }
        class.prob <- nd / rs
      } else {
        stop("'newdata' for FCGDINA must be an M x D binary alpha matrix ",
             "or an M x 2^D class-probability matrix.", call. = FALSE)
      }
    }
    return(class.prob %*% prob.class)
  }

  theta <- if (is.null(newdata)) as.matrix(object$theta$est) else as.matrix(newdata)
  switch(cls,
    MIRT   = model.MIRT(theta, as.matrix(object$par$est)),
    MGPCM  = model.MGPCM(theta, as.matrix(object$par$est)),
    MGGUM  = model.MGGUM(theta, as.matrix(object$par$est)),
    FCMIRT = model.FCMIRT(theta, as.matrix(object$par$est),
                          object$patterns.total, object$patterns),
    FCDCM  = {
      alpha_prob <- fcdcm_alpha_profile_prob(
        theta = theta,
        delta1 = object$delta$est[, "delta1"],
        delta0 = object$delta$est[, "delta0"],
        alpha.patterns = object$alpha.patterns
      )
      prob_per_pattern <- model.FCDCM(
        zeta = object$zeta.patterns,
        par = object$par$est,
        patterns = object$patterns
      )
      alpha_prob %*% prob_per_pattern
    },
    FCGGUM = model.FCGGUM(theta, as.matrix(object$par$est),
                          object$patterns.total, object$patterns),
    TIRT   = {
      par <- as.matrix(object$par$est)
      D <- ncol(par) - 1L
      lambda <- unlist(lapply(seq_len(nrow(par)), function(x) {
        d <- which(object$Q.matrix[x, ] != 0)[1]
        return(par[x, d])
      }))
      psi2 <- par[, D + 1L]
      model.TIRT(theta, lambda, psi2, object$gamma.matrix$est,
                 object$Q.matrix, object$pairs.matrix)
    }
  )
}

.fc_predict_response_binary <- function(prob) {
  matrix(as.integer(runif(length(prob)) < prob), nrow = nrow(prob),
         ncol = ncol(prob), dimnames = dimnames(prob))
}

.fc_predict_response_polytomous <- function(prob, length.poly) {
  N <- nrow(prob)
  I <- length(length.poly)
  resp <- matrix(0L, N, I)
  idx <- 1L
  for (i in seq_len(I)) {
    Ki <- length.poly[i]
    prob_i <- prob[, idx:(idx + Ki - 1L), drop = FALSE]
    cum_prob <- t(apply(prob_i, 1L, cumsum))
    u <- runif(N)
    resp[, i] <- rowSums(u > cum_prob)
    idx <- idx + Ki
  }
  resp
}

.fc_predict_response_fc <- function(prob, patterns) {
  N <- nrow(prob)
  B <- length(patterns)
  n_cat <- vapply(patterns, nrow, integer(1L))
  resp <- matrix(0L, N, B)
  idx <- 1L
  for (b in seq_len(B)) {
    Kb <- n_cat[b]
    prob_b <- prob[, idx:(idx + Kb - 1L), drop = FALSE]
    cum_prob <- t(apply(prob_b, 1L, cumsum))
    u <- runif(N)
    resp[, b] <- rowSums(u > cum_prob) + 1L
    idx <- idx + Kb
  }
  resp
}


# ===========================================================================
# MIRT
# ===========================================================================

#' @describeIn predict Predict response probabilities or simulated
#'   responses from \code{MIRT} objects.
#' @method predict MIRT
#' @export
predict.MIRT <- function(object, newdata = NULL,
                         type = c("probability", "response"), ...) {
  type <- match.arg(type)
  prob <- .fc_predict_prob(object, newdata)
  switch(type,
    probability = prob,
    response    = .fc_predict_response_binary(prob)
  )
}

# ===========================================================================
# MGPCM
# ===========================================================================

#' @describeIn predict Predict response probabilities from \code{MGPCM} objects.
#' @method predict MGPCM
#' @export
predict.MGPCM <- function(object, newdata = NULL,
                          type = c("probability", "response"), ...) {
  type <- match.arg(type)
  prob <- .fc_predict_prob(object, newdata)
  switch(type,
    probability = prob,
    response    = .fc_predict_response_polytomous(prob, object$length.poly)
  )
}

# ===========================================================================
# MGGUM
# ===========================================================================

#' @describeIn predict Predict response probabilities from \code{MGGUM} objects.
#' @method predict MGGUM
#' @export
predict.MGGUM <- function(object, newdata = NULL,
                          type = c("probability", "response"), ...) {
  type <- match.arg(type)
  prob <- .fc_predict_prob(object, newdata)
  switch(type,
    probability = prob,
    response    = .fc_predict_response_polytomous(prob, object$length.poly)
  )
}

# ===========================================================================
# FCMIRT
# ===========================================================================

#' @describeIn predict Predict response probabilities from \code{FCMIRT} objects.
#' @method predict FCMIRT
#' @export
predict.FCMIRT <- function(object, newdata = NULL,
                           type = c("probability", "response"), ...) {
  type <- match.arg(type)
  prob <- .fc_predict_prob(object, newdata)
  switch(type,
    probability = prob,
    response    = .fc_predict_response_fc(prob, object$patterns)
  )
}

# ===========================================================================
# FCDCM
# ===========================================================================

#' @describeIn predict Predict response probabilities from \code{FCDCM} objects.
#' @method predict FCDCM
#' @export
predict.FCDCM <- function(object, newdata = NULL,
                          type = c("probability", "response"), ...) {
  type <- match.arg(type)
  prob <- .fc_predict_prob(object, newdata)
  switch(type,
    probability = prob,
    response    = .fc_predict_response_binary(prob)
  )
}

# ===========================================================================
# FCGDINA
# ===========================================================================

#' @describeIn predict Predict response probabilities from \code{FCGDINA} objects.
#'   For \code{newdata}, provide either binary alpha profiles or class
#'   probability rows.
#' @method predict FCGDINA
#' @export
predict.FCGDINA <- function(object, newdata = NULL,
                            type = c("probability", "response"), ...) {
  type <- match.arg(type)
  prob <- .fc_predict_prob(object, newdata)
  switch(type,
    probability = prob,
    response    = .fc_predict_response_fc(prob, object$patterns)
  )
}

# ===========================================================================
# FCGGUM
# ===========================================================================

#' @describeIn predict Predict response probabilities from \code{FCGGUM} objects.
#' @method predict FCGGUM
#' @export
predict.FCGGUM <- function(object, newdata = NULL,
                           type = c("probability", "response"), ...) {
  type <- match.arg(type)
  prob <- .fc_predict_prob(object, newdata)
  switch(type,
    probability = prob,
    response    = .fc_predict_response_fc(prob, object$patterns)
  )
}

# ===========================================================================
# TIRT
# ===========================================================================

#' @describeIn predict Predict response probabilities from \code{TIRT} objects.
#' @method predict TIRT
#' @export
predict.TIRT <- function(object, newdata = NULL,
                         type = c("probability", "response"), ...) {
  type <- match.arg(type)
  prob <- .fc_predict_prob(object, newdata)
  switch(type,
    probability = prob,
    response    = .fc_predict_response_binary(prob)
  )
}
