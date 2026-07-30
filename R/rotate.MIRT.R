# ===========================================================================
# rotate.MIRT
# ===========================================================================

#' @describeIn rotate Applies rotation to the posterior mean of loadings
#'   only from \code{MIRT} objects. Loadings (a) and factor correlation
#'   matrix (Corr) are updated to preserve model-implied covariances.
#'   For GPArotation methods, theta is re-standardized to a standard
#'   normal scale after rotation; for Promax, no post-rotation
#'   standardization is applied.
#' @method rotate MIRT
#' @export
rotate.MIRT <- function(object,
                        method,
                        vis          = TRUE,
                        target.a     = NULL,
                        target.theta = NULL,
                        pst.W        = NULL,
                        lp.p         = 1,
                        lp.gpaiter   = 5,
                        promax_m     = 4,
                        ...) {

  call <- match.call()
  user_args <- list(...)

  if (!is.character(method) || length(method) != 1L ||
      is.na(method) || !nzchar(method)) {
    stop("'method' must be a non-empty character scalar.", call. = FALSE)
  }

  D <- ncol(object$par$est) - 3L
  I <- nrow(object$par$est)

  rotation_arguments <- c(
    list(
      object = object,
      method = method,
      vis = vis,
      target.a = target.a,
      target.theta = target.theta,
      pst.W = pst.W,
      lp.p = lp.p,
      lp.gpaiter = lp.gpaiter,
      promax_m = promax_m
    ),
    user_args
  )
  new_rotation_result <- function(out) {
    structure(
      list(object = out, call = call, arguments = rotation_arguments),
      class = "RotateMIRT"
    )
  }

  if (D == 1L) {
    if (vis) message("D = 1: rotation is not applicable. Returning original fit object.")
    object$rotated         <- FALSE
    object$rotation_method <- method
    return(new_rotation_result(object))
  }

  model <- tolower(object$arguments$model)
  if (identical(model, "m1pl")) {
    stop("Rotation is not defined for multidimensional M1PL models with fixed loadings.",
         call. = FALSE)
  }

  Q.matrix <- object$Q.matrix
  if (!is.null(Q.matrix) && any(Q.matrix == 0)) {
    stop("Rotation is not allowed when Q.matrix contains zeros (confirmatory MIRT). ",
         "The Q.matrix specifies a confirmatory structure that rotation would destroy.",
         call. = FALSE)
  }

  use_target_a     <- !is.null(target.a)
  use_target_theta <- !is.null(target.theta)
  is_target_method <- method %in% c("targetT", "targetQ", "pstT", "pstQ")

  if (is_target_method) {
    if (use_target_a && use_target_theta) {
      if (vis) warning("Both 'target.a' and 'target.theta' provided. 'target.a' takes precedence; 'target.theta' will be ignored.")
      use_target_theta <- FALSE
    }
    if (!use_target_a && !use_target_theta) {
      stop(sprintf("Method '%s' requires either 'target.a' or 'target.theta'.", method),
           call. = FALSE)
    }
    if (use_target_theta && method %in% c("pstT", "pstQ")) {
      stop(sprintf("Method '%s' requires 'target.a'. 'target.theta' is not supported for PST.", method),
           call. = FALSE)
    }
    if (method %in% c("pstT", "pstQ") && is.null(pst.W)) {
      stop(sprintf("Method '%s' requires a weight matrix 'pst.W'.", method),
           call. = FALSE)
    }
    if (use_target_a && ncol(target.a) != D) {
      stop(sprintf("target.a has %d columns but model D = %d.", ncol(target.a), D),
           call. = FALSE)
    }
    if (use_target_theta && ncol(target.theta) != D) {
      stop(sprintf("target.theta has %d columns but model D = %d.", ncol(target.theta), D),
           call. = FALSE)
    }
  } else {
    if (use_target_a || use_target_theta) {
      if (vis) warning(sprintf("target.a/target.theta provided but method '%s' does not use targets. Ignoring.", method))
      use_target_a <- FALSE
      use_target_theta <- FALSE
    }
  }

  A_mean     <- object$par$est[, 1:D, drop = FALSE]
  theta_mean <- object$theta$est
  Corr_mean  <- object$Corr$est

  if (tolower(method) == "promax") {
    if (use_target_a || use_target_theta) {
      if (vis) warning("target.a/target.theta provided but method 'promax' does not use targets. Ignoring.")
    }

    if (vis) message(sprintf("Rotating posterior mean using 'promax' (m = %d)...", promax_m))

    rot_result_raw <- tryCatch(
      stats::promax(A_mean, m = promax_m),
      error = function(e) {
        stop(sprintf("Rotation of posterior mean failed: %s", conditionMessage(e)),
             call. = FALSE)
      }
    )

    T_mat <- rot_result_raw$rotmat
    standardize <- FALSE

  } else {
    base_args <- list(
      Tmat         = diag(D),
      normalize    = FALSE,
      eps          = 1e-5,
      maxit        = 1000L,
      randomStarts = 0L
    )

    method_specific <- switch(method,
                              "targetT" = , "targetQ" = list(Target = NULL, L = NULL),
                              "pstT"    = , "pstQ"    = list(W = pst.W, Target = NULL, L = NULL),
                              "oblimin"               = list(gam = 0),
                              "quartimin"             = list(),
                              "Varimax"               = list(),
                              "quartimax"             = list(),
                              "simplimax"             = list(k = I),
                              "geominT" = , "geominQ" = ,
                              "bigeominT" = , "bigeominQ" = list(delta = 0.01),
                              "cfT" = , "cfQ"         = list(kappa = 0),
                              "equamax"               = list(kappa = D / (2 * I)),
                              "parsimax"              = list(kappa = (D - 1) / (D + I - 2)),
                              "lpT" = , "lpQ"         = list(p = lp.p, gpaiter = lp.gpaiter),
                              "entropy" = , "oblimax" = ,
                              "bentlerT" = , "bentlerQ" = ,
                              "tandemI" = , "tandemII" = ,
                              "infomaxT" = , "infomaxQ" = ,
                              "mccammon" = , "varimin" = ,
                              "bifactorT" = , "bifactorQ" = list(),
                              list()
    )

    call_template <- utils::modifyList(base_args, method_specific)
    call_template <- utils::modifyList(call_template, user_args)
    call_template <- call_template[!sapply(call_template, is.null)]

    rot_fn <- tryCatch(
      get(method, mode = "function", envir = asNamespace("GPArotation")),
      error = function(e) {
        stop(sprintf("Rotation method '%s' not found in GPArotation package.", method),
             call. = FALSE)
      }
    )

    iter_call <- call_template

    if (is_target_method) {
      if (use_target_a) {
        iter_call$Target <- target.a
      } else if (use_target_theta) {
        svd_res <- svd(crossprod(target.theta, theta_mean))
        bridge  <- svd_res$v %*% t(svd_res$u)
        iter_call$Target <- A_mean %*% bridge
      }
    }

    iter_call$A <- A_mean

    if (vis) message(sprintf("Rotating posterior mean using '%s'...", method))

    rot_result <- tryCatch(
      do.call(rot_fn, iter_call),
      error = function(e) {
        stop(sprintf("Rotation of posterior mean failed: %s", conditionMessage(e)),
             call. = FALSE)
      }
    )

    if (is.null(rot_result) || is.null(rot_result$Th)) {
      stop("Rotation of posterior mean failed: no transformation matrix returned.",
           call. = FALSE)
    }

    T_mat <- if (isTRUE(rot_result$orthogonal)) {
      rot_result$Th
    } else {
      solve(t(rot_result$Th))
    }
    standardize <- TRUE
  }

  rotated <- .fc_transform_rotation(
    a = A_mean,
    theta = theta_mean,
    Corr = Corr_mean,
    Tmat = T_mat,
    standardize = standardize
  )

  par   <- object$par$est
  par[, seq_len(D)] <- rotated$a
  dimnames(par) <- dimnames(object$par$est)

  out <- object
  out$par$est   <- par
  out$theta$est <- rotated$theta
  out$Corr$est  <- rotated$Corr

  out$rotated         <- TRUE
  out$rotation_method <- method

  new_rotation_result(out)
}
