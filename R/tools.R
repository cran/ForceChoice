# ---- Control helpers ----------------------------------------------------------

get_ctrl <- function(name, default, control) {
  if (!is.null(control[[name]])) control[[name]] else default
}


# ---- Array helpers ------------------------------------------------------------

ensure_3d <- function(x) {
  if (is.null(x)) return(NULL)
  nd <- length(dim(x))
  if (nd == 0L) {
    # 1D vector -> (iter, 1, 1)
    dim(x) <- c(length(x), 1L, 1L)
  } else if (nd == 2L) {
    # 2D matrix  -> (iter, rows, 1)
    dim(x) <- c(dim(x), 1L)
  }
  # nd == 3L: already correct, return as-is
  return(x)
}


# ---- FC type normalisation ----------------------------------------------------

#' Normalise and recycle forced-choice type vector
#'
#' @param fc.type Character vector (or \code{NULL}) describing the FC format of
#'   each block.  Recycled to length \code{N.block} and uppercased.
#' @param N.block Number of blocks.
#' @param default Default FC type used when \code{fc.type} is \code{NULL}.
#' @return A character vector of length \code{N.block}.
#' @keywords internal
normalize_fc_type <- function(fc.type, N.block, default = "RANK") {
  N.block <- check_integer_scalar(N.block, "N.block")
  if (is.null(fc.type)) fc.type <- default
  if (length(fc.type) == 1L) {
    fc.type <- rep(fc.type, N.block)
  } else if (length(fc.type) != N.block) {
    stop("'fc.type' must have length 1 or the number of blocks (",
         N.block, ").", call. = FALSE)
  }
  fc.type <- toupper(trimws(as.character(fc.type)))
  valid <- c("RANK", "MOLE", "PICK")
  if (anyNA(fc.type) || any(!fc.type %in% valid)) {
    stop("'fc.type' must contain only RANK, MOLE, or PICK.",
         call. = FALSE)
  }
  fc.type
}


# ---- Permutations  ------------------------------------------------------------

get_permutations <- function(vector, r) {
  n <- length(vector)
  if (r > n) stop("Selected element count 'r' cannot exceed the vector length.", call. = FALSE)
  vector <- as.integer(vector)
  if (anyNA(vector)) stop("'vector' must not contain NA values.", call. = FALSE)
  cpp_get_permutations(vector, as.integer(r))
}


# ---- Correlation matrix validation -------------------------------------------

#' Validate a user-supplied correlation matrix
#'
#' Checks dimensions, symmetry, unit diagonal, and positive semi-definiteness.
#' Returns a valid correlation matrix (or \code{diag(D)} when \code{NULL}).
#'
#' @param Corr A D x D matrix or \code{NULL}.
#' @param D Expected number of dimensions.
#' @param label Character label for error messages.
#' @return A valid D x D correlation matrix.
#' @keywords internal
validate_corr_matrix <- function(Corr, D, label = "'Corr'") {
  if (is.null(Corr)) return(diag(D))
  if (!is.matrix(Corr) || nrow(Corr) != D || ncol(Corr) != D)
    stop(label, " must be a ", D, "x", D, " matrix", call. = FALSE)
  if (!isSymmetric(unname(Corr), tol = 1e-8))
    stop(label, " must be symmetric", call. = FALSE)
  if (any(abs(diag(Corr) - 1) > 1e-8))
    stop(label, " must be a correlation matrix with diagonal 1", call. = FALSE)
  eig <- eigen(Corr, symmetric = TRUE, only.values = TRUE)$values
  if (any(eig < -1e-10))
    stop(label, " is not positive semi-definite", call. = FALSE)
  Corr
}


# ---- Rhat extraction ----------------------------------------------------------

#' Extract R-hat values for a matrix parameter from a Stan summary
#'
#' @param stan.sum Stan summary matrix (rownames encode parameter indices).
#' @param param_name Parameter name prefix (e.g. \code{"theta"}).
#' @param nrow,ncol Dimensions of the output matrix.
#' @return A numeric matrix of R-hat values.
#' @keywords internal
extract_rhat_matrix <- function(stan.sum, param_name, nrow, ncol) {
  cpp_extract_rhat_matrix(
    rhat       = stan.sum[, "Rhat"],
    row_names  = rownames(stan.sum),
    param_name = param_name,
    nrow_out   = as.integer(nrow),
    ncol_out   = as.integer(ncol)
  )
}

#' Extract R-hat values for a vector parameter from a Stan summary
#'
#' @param stan.sum Stan summary matrix.
#' @param param_name Parameter name prefix (e.g. \code{"delta1"}).
#' @param n Expected length of the output vector.
#' @return A numeric vector of R-hat values (NA where not found).
#' @keywords internal
extract_rhat_vector <- function(stan.sum, param_name, n) {
  cpp_extract_rhat_vector(
    rhat       = stan.sum[, "Rhat"],
    row_names  = rownames(stan.sum),
    param_name = param_name,
    n_out      = as.integer(n)
  )
}


# ---- Parameter-counting helpers -----------------------------------------------

#' Degrees of freedom for a correlation matrix
#'
#' @param D Number of dimensions.
#' @param include_corr Logical; if \code{FALSE} returns 0.
#' @return Integer, the number of free parameters in a D x D correlation matrix.
#' @keywords internal
free_corr_npar <- function(D, include_corr = TRUE) {
  if (isTRUE(include_corr) && D > 1L) D * (D - 1L) / 2L else 0L
}

#' Free-parameter mask for MIRT models
#'
#' @param model_type Integer 2, 3, or 4.
#' @param Q.matrix I x D Q-matrix.
#' @return An I x (D+3) logical matrix.
#' @keywords internal
mirt_par_free_mask <- function(model_type, Q.matrix) {
  I <- nrow(Q.matrix)
  D <- ncol(Q.matrix)
  free <- matrix(FALSE, I, D + 3L)
  if (model_type >= 2L) {
    free[, seq_len(D)] <- Q.matrix == 1
  }
  free[, D + 1L] <- TRUE                                # b always free
  if (model_type >= 3L) {
    free[, D + 2L] <- TRUE                              # c
  }
  if (model_type == 4L) {
    free[, D + 3L] <- TRUE                              # d
  }
  colnames(free) <- c(paste0("a", seq_len(D)), "b", "c", "d")
  free
}

#' Free-parameter mask for MGPCM models
#'
#' @param Q.matrix I x D Q-matrix.
#' @param length.poly Integer vector of category counts per item.
#' @param max_poly Maximum number of categories (default \code{max(length.poly)}).
#' @return An I x (D + max_poly) logical matrix.
#' @keywords internal
mgpcm_par_free_mask <- function(Q.matrix, length.poly,
                                max_poly = max(length.poly)) {
  I <- nrow(Q.matrix)
  D <- ncol(Q.matrix)
  free <- matrix(FALSE, I, D + max_poly)
  free[, seq_len(D)] <- Q.matrix == 1
  for (i in seq_len(I)) {
    Ki <- length.poly[i]
    if (Ki > 1L) {
      free[i, D + 2:Ki] <- TRUE
    }
  }
  colnames(free) <- c(paste0("a", seq_len(D)),
                      paste0("d", 0:(max_poly - 1L)))
  free
}

#' Free-parameter mask for MGGUM models
#'
#' @param Q.matrix I x D Q-matrix.
#' @param length.poly Integer vector of category counts per item.
#' @param max_poly Maximum number of categories.
#' @return An I x (2D + max_poly) logical matrix.
#' @keywords internal
mggum_par_free_mask <- function(Q.matrix, length.poly,
                                max_poly = max(length.poly)) {
  I <- nrow(Q.matrix)
  D <- ncol(Q.matrix)
  free <- matrix(FALSE, I, D + D + max_poly)
  active <- Q.matrix != 0
  free[, seq_len(D)]         <- active
  free[, D + seq_len(D)]     <- active
  for (i in seq_len(I)) {
    Ki <- length.poly[i]
    if (Ki > 1L) {
      free[i, D + D + 2:Ki] <- TRUE
    }
  }
  colnames(free) <- c(
    paste0("a",     seq_len(D)),
    paste0("delta", seq_len(D)),
    paste0("tau",   0:(max_poly - 1L))
  )
  free
}

#' Free-parameter mask for TIRT models
#'
#' @param Q.matrix I x D Q-matrix.
#' @param lambda.free Logical vector of length I.
#' @param psi.free Logical vector of length I.
#' @return An I x (D+1) logical matrix.
#' @keywords internal
tirt_par_free_mask <- function(Q.matrix, lambda.free, psi.free) {
  I <- nrow(Q.matrix)
  D <- ncol(Q.matrix)
  free <- matrix(FALSE, I, D + 1L)
  for (i in seq_len(I)) {
    if (isTRUE(lambda.free[i])) {
      active <- which(Q.matrix[i, ] != 0)
      if (length(active) > 0L) {
        free[i, active[1L]] <- TRUE
      }
    }
  }
  free[, D + 1L] <- psi.free
  colnames(free) <- c(paste0("Dim.", seq_len(D)), "psi2")
  free
}
