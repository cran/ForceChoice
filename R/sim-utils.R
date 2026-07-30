# ---- Shared utilities for sim.data.* functions ---------------------------------
#
# This file contains helper functions that are used across multiple
# sim.data.* data-simulation functions.  They are tagged @keywords internal
# and are not exported.


# ---- Integer scalar validation ------------------------------------------------

#' Validate a scalar integer parameter
#'
#' Checks that \code{x} is a single finite value coercible to integer and
#' at least \code{min.val}.  Used by forced-choice \code{sim.data.*}
#' functions to validate \code{N.person}, \code{N.block}, \code{I.block},
#' and \code{D}.
#'
#' @param x Value to check.
#' @param name Parameter name (used in the error message).
#' @param min.val Minimum allowed value (default 1L).
#' @return \code{as.integer(x)} on success; raises an error otherwise.
#' @keywords internal
check_integer_scalar <- function(x, name, min.val = 1L) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) ||
      x < min.val || x != floor(x)) {
    stop("'", name, "' must be an integer >= ", min.val, ".",
         call. = FALSE)
  }
  as.integer(x)
}

normalize_length_poly <- function(length.poly, I) {
  if (length(length.poly) == 1L) {
    length.poly <- rep(length.poly, I)
  } else if (length(length.poly) != I) {
    stop("'length.poly' must have length 1 or 'I'.", call. = FALSE)
  }
  if (!is.numeric(length.poly) || anyNA(length.poly) ||
      any(!is.finite(length.poly)) || any(length.poly != floor(length.poly)) ||
      any(length.poly < 2L)) {
    stop("'length.poly' must contain integers >= 2.", call. = FALSE)
  }
  as.integer(length.poly)
}


# ---- Theta generation ---------------------------------------------------------

#' Generate latent trait vectors from a multivariate normal distribution
#'
#' When \code{D > 1}, uses \code{MASS::mvrnorm} with the supplied correlation
#' matrix.  When \code{D == 1}, uses standard normal draws.
#'
#' @param N Number of persons.
#' @param D Number of latent dimensions.
#' @param Corr A \eqn{D \times D} correlation matrix.
#' @return An \eqn{N \times D} numeric matrix of latent trait values.
#' @keywords internal
generate_theta_mvn <- function(N, D, Corr) {
  if (D > 1L) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for multivariate normal theta generation.", call. = FALSE)
    }
    MASS::mvrnorm(n = N, mu = rep(0, D), Sigma = Corr)
  } else {
    matrix(rnorm(N, 0, 1), ncol = 1L)
  }
}


# ---- Factor rotation ----------------------------------------------------------

# Apply a loading-space transformation consistently to loadings, row-wise
# latent scores, and their column covariance matrix.  If A* = A T, preserving
# theta %*% t(A) requires theta* = theta T^{-T}; consequently the latent
# covariance becomes T^{-1} Corr T^{-T}.
#
# @param a I x D loading matrix.
# @param theta N x D matrix whose rows are latent scores.
# @param Corr D x D latent covariance/correlation matrix.
# @param Tmat D x D loading transformation matrix.
# @param standardize Whether to rescale the transformed latent dimensions to
#   unit variance while applying the inverse rescaling to theta and a.
# @return A list containing transformed a, theta, Corr, and the effective Tmat.
# @noRd
.fc_transform_rotation <- function(a, theta, Corr, Tmat,
                                   standardize = FALSE) {
  a <- as.matrix(a)
  theta <- as.matrix(theta)
  Corr <- as.matrix(Corr)
  Tmat <- as.matrix(Tmat)

  D <- ncol(a)
  if (ncol(theta) != D ||
      !identical(dim(Corr), c(D, D)) ||
      !identical(dim(Tmat), c(D, D))) {
    stop("'a', 'theta', 'Corr', and 'Tmat' have incompatible dimensions.",
         call. = FALSE)
  }
  if (anyNA(a) || anyNA(theta) || anyNA(Corr) || anyNA(Tmat) ||
      any(!is.finite(a)) || any(!is.finite(theta)) ||
      any(!is.finite(Corr)) || any(!is.finite(Tmat))) {
    stop("Rotation inputs must contain only finite numeric values.",
         call. = FALSE)
  }

  T_inv <- tryCatch(
    solve(Tmat),
    error = function(e) {
      stop("Rotation transformation matrix is singular: ",
           conditionMessage(e), call. = FALSE)
    }
  )

  a_rot <- a %*% Tmat
  theta_rot <- theta %*% t(T_inv)
  Corr_rot <- T_inv %*% Corr %*% t(T_inv)
  Corr_rot <- (Corr_rot + t(Corr_rot)) / 2
  T_effective <- Tmat

  if (isTRUE(standardize)) {
    variances <- diag(Corr_rot)
    if (any(!is.finite(variances)) ||
        any(variances <= .Machine$double.eps)) {
      stop("Rotated latent variances must be finite and positive.",
           call. = FALSE)
    }

    scale_mat <- diag(sqrt(variances), nrow = D)
    scale_inv <- diag(1 / sqrt(variances), nrow = D)
    a_rot <- a_rot %*% scale_mat
    theta_rot <- theta_rot %*% scale_inv
    Corr_rot <- scale_inv %*% Corr_rot %*% scale_inv
    Corr_rot <- (Corr_rot + t(Corr_rot)) / 2
    diag(Corr_rot) <- 1
    T_effective <- Tmat %*% scale_mat
  }

  dimnames(a_rot) <- dimnames(a)
  dimnames(theta_rot) <- dimnames(theta)
  dimnames(Corr_rot) <- dimnames(Corr)
  dimnames(T_effective) <- dimnames(Tmat)

  list(
    a = a_rot,
    theta = theta_rot,
    Corr = Corr_rot,
    Tmat = T_effective
  )
}

#' Apply factor rotation to discrimination matrix and latent traits
#'
#' Supports \code{"promax"} (via \code{stats::promax}) and GPArotation
#' methods (\code{"varimax"}, \code{"oblimin"}, etc.).  Transforms the
#' discrimination matrix \code{a}, the latent trait matrix \code{theta},
#' and the correlation matrix \code{Corr} accordingly.
#'
#' @param a \eqn{I \times D} discrimination (loading) matrix.
#' @param theta \eqn{N \times D} latent trait matrix.
#' @param Corr \eqn{D \times D} correlation matrix.
#' @param rotate Rotation method name (character).
#' @param D Number of dimensions.
#' @param promax_m Power parameter for \code{"promax"} (default 4).
#' @return A named list with components \code{a}, \code{theta}, \code{Corr}.
#' @keywords internal
apply_rotation_promax <- function(a, theta, Corr, rotate, D, promax_m = 4) {
  if (is.null(rotate) || D <= 1L) {
    return(list(a = a, theta = theta, Corr = Corr))
  }

  if (tolower(rotate) == "promax") {
    if (!requireNamespace("stats", quietly = TRUE)) {
      stop("Package 'stats' is required for promax rotation.", call. = FALSE)
    }
    rot.res <- stats::promax(a, m = promax_m)
    T.mat <- rot.res$rotmat
  } else {
    if (rotate == "Varimax") {
      rotate <- "varimax"
    }
    orth.methods <- c("varimax", "quartimax")
    obl.methods  <- c("oblimin", "quartimin", "oblimax",
                      "entropy", "simplimax", "geominQ")
    if (rotate %in% orth.methods) {
      if (!requireNamespace("GPArotation", quietly = TRUE)) {
        stop("Package 'GPArotation' is required for rotation method '", rotate, "'.", call. = FALSE)
      }
      rot.res <- GPArotation::GPForth(a, method = rotate)
      T.mat <- rot.res$Th
    } else if (rotate %in% obl.methods) {
      if (!requireNamespace("GPArotation", quietly = TRUE)) {
        stop("Package 'GPArotation' is required for rotation method '", rotate, "'.", call. = FALSE)
      }
      rot.res <- GPArotation::GPFoblq(a, method = rotate)
      T.mat <- solve(t(rot.res$Th))
    } else {
      warning("Unknown rotation method '", rotate, "'. Skipping rotation.")
      return(list(a = a, theta = theta, Corr = Corr))
    }
  }

  rotated <- .fc_transform_rotation(
    a = a,
    theta = theta,
    Corr = Corr,
    Tmat = T.mat,
    standardize = TRUE
  )
  rotated[c("a", "theta", "Corr")]
}


# ---- Forced-choice block / Q-matrix resolution --------------------------------

#' Resolve forced-choice block and Q-matrix configuration
#'
#' Shared validation and auto-generation logic used by
#' \code{sim.data.FCMIRT}, \code{sim.data.FCGGUM}, \code{sim.data.FCGDINA},
#' and \code{sim.data.TIRT}.
#'
#' @param Q.matrix User-supplied Q-matrix or \code{NULL}.
#' @param block.items User-supplied block-item list or \code{NULL}.
#' @param I.block Default items-per-block (scalar).
#' @param N.block Default number of blocks (scalar).
#' @param D Default number of dimensions (scalar).
#' @param q.valid.values Allowed Q-matrix values (e.g. \code{c(0, 1)} for
#'   MIRT/DCM, \code{c(-1, 0, 1)} for GGUM/TIRT).
#' @param q.row.check Function applied to each row of Q.matrix; must return
#'   \code{TRUE} for every row.  Example: \code{function(r) any(r != 0)} or
#'   \code{function(r) sum(abs(r)) == 1L}.
#'
#' @return A named list with components \code{I.states}, \code{N.block},
#'   \code{I.block}, \code{block.items}, \code{Q.matrix}, \code{D}.
#'   \code{Q.matrix} is \code{NULL} if it was not supplied and should be
#'   generated later (e.g. by \code{sim.Q.matrix.FC} or \code{sim.data.Q.CDM}).
#' @keywords internal
resolve_fc_blocks <- function(Q.matrix, block.items, I.block, N.block, D,
                              q.valid.values, q.row.check) {
  if (!is.null(Q.matrix)) {
    Q.matrix <- as.matrix(Q.matrix)
    if (anyNA(Q.matrix) || !all(Q.matrix %in% q.valid.values)) {
      stop("'Q.matrix' must contain only ",
           paste(q.valid.values, collapse = ", "), ".", call. = FALSE)
    }
    if (!all(apply(Q.matrix, 1L, q.row.check))) {
      stop("Each row of 'Q.matrix' violates the required row constraint.", call. = FALSE)
    }
    storage.mode(Q.matrix) <- "numeric"
    I.states <- nrow(Q.matrix)
    D        <- ncol(Q.matrix)
  } else {
    I.states <- N.block * I.block
  }

  if (!is.null(block.items)) {
    if (!is.list(block.items)) {
      stop("'block.items' must be a list.", call. = FALSE)
    }
    N.block <- length(block.items)
    if (N.block < 1L) {
      stop("'block.items' must contain at least one block.", call. = FALSE)
    }

    block.lengths <- vapply(block.items, length, integer(1L))
    if (any(block.lengths < 2L)) {
      stop("Each block in 'block.items' must contain at least two items.", call. = FALSE)
    }
    I.block <- if (length(unique(block.lengths)) == 1L) {
      block.lengths[1L]
    } else {
      block.lengths
    }

    all.items <- unlist(block.items, use.names = FALSE)
    if (anyDuplicated(all.items)) {
      stop("'block.items' contains duplicate items.", call. = FALSE)
    }

    if (is.null(Q.matrix)) {
      I.states <- length(all.items)
    }
    if (!setequal(all.items, seq_len(I.states))) {
      stop("'block.items' must contain exactly all item indices from 1 to the number of items.", call. = FALSE)
    }
  } else {
    if (is.null(Q.matrix)) {
      I.states <- N.block * I.block
    } else {
      if (I.states %% I.block != 0L) {
        stop("Number of items in 'Q.matrix' is not divisible by I.block.", call. = FALSE)
      }
      N.block <- as.integer(I.states / I.block)
    }
    block.items <- split(seq_len(I.states),
                         rep(seq_len(N.block), each = I.block))
  }

  list(I.states   = I.states,
       N.block    = N.block,
       I.block    = I.block,
       block.items = block.items,
       Q.matrix   = Q.matrix,
       D          = D)
}


# ---- Forced-choice permutation patterns ---------------------------------------

#' Generate forced-choice permutation patterns
#'
#' Produces the full-ranking patterns (\code{patterns.total}) and the
#' observed patterns (\code{patterns}) for each block, based on the
#' forced-choice response type.
#'
#' @param block.items List of integer vectors (item indices per block).
#' @param fc.type Character vector of length \code{length(block.items)}
#'   with values \code{"RANK"}, \code{"MOLE"}, or \code{"PICK"}.
#' @return A named list with components \code{patterns.total} (list of
#'   full-ranking matrices) and \code{patterns} (list of observed-pattern
#'   matrices).
#' @keywords internal
generate_fc_permutation_patterns <- function(block.items, fc.type) {
  N.block <- length(block.items)

  patterns.total <- lapply(block.items, function(x) {
    get_permutations(x, length(x))
  })

  patterns <- vector("list", N.block)
  for (bb in seq_len(N.block)) {
    patterns[[bb]] <- get_permutations(
      block.items[[bb]],
      switch(fc.type[bb],
             RANK = length(block.items[[bb]]),
             MOLE = 2L,
             PICK = 1L)
    )
  }

  list(patterns.total = patterns.total, patterns = patterns)
}


# ---- Forced-choice response sampling ------------------------------------------

#' Sample forced-choice responses from block-level pattern probabilities
#'
#' Uses inverse-CDF multinomial sampling (identical to the C++ implementation
#' in \code{forced_choice_from_agree}) to draw a discrete pattern index for
#' each person-block combination.
#'
#' @param prob \eqn{N \times \sum_B P_b} matrix of block-level pattern
#'   probabilities (stacked by block).
#' @param patterns List of observed-pattern matrices (one per block).
#' @param N.person Number of persons.
#' @param N.block Number of blocks.
#' @return A named list with \code{response} (\eqn{N \times B} integer
#'   matrix of pattern indices) and \code{data} (\eqn{N \times B} character
#'   matrix of ranking strings like \code{"1>3>2"}).
#' @keywords internal
sample_fc_responses <- function(prob, patterns, N.person, N.block) {
  response <- matrix(NA_integer_, N.person, N.block)
  data     <- matrix(NA_character_, N.person, N.block)

  idx <- 0L
  for (bb in seq_len(N.block)) {
    patterns.cur <- patterns[[bb]]
    bR <- nrow(patterns.cur)
    cols <- idx + seq_len(bR)
    for (pp in seq_len(N.person)) {
      prob.bp <- cumsum(prob[pp, cols])
      y <- which(runif(1, 0, 1) <= prob.bp)[1L]
      if (is.na(y)) {
        y <- bR
      }
      response[pp, bb] <- y
      data[pp, bb] <- paste(patterns.cur[y, ], collapse = ">")
    }
    idx <- idx + bR
  }

  list(response = response, data = data)
}


# ---- Pattern key conversion ---------------------------------------------------

#' Convert binary matrix rows to string keys
#'
#' Collapses each row of a binary matrix to a single string by pasting
#' the entries without separator (e.g. \code{c(0,1,0)} becomes \code{"010"}).
#'
#' @param x A numeric or integer matrix.
#' @return A character vector of length \code{nrow(x)}.
#' @keywords internal
pattern_key <- function(x) {
  x <- as.matrix(x)
  apply(x, 1L, paste0, collapse = "")
}


# ---- Probability clamping -----------------------------------------------------

#' Clamp probabilities away from 0 and 1
#'
#' Ensures all values lie in \eqn{[eps, 1-eps]} to avoid \code{log(0)}
#' or division-by-zero issues in downstream computations.
#'
#' @param p Numeric vector or matrix of probabilities.
#' @param eps Small positive value (default \code{1e-12}).
#' @return Numeric vector or matrix with the same dimensions as \code{p}.
#' @keywords internal
clip_prob <- function(p, eps = 1e-12) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}


# ---- CDM item response probabilities ------------------------------------------

#' Build CDM item response probabilities for attribute mastery patterns
#'
#' Computes the probability of endorsing an item for each attribute mastery
#' pattern under the specified CDM condensation rule (DINA, DINO, ACDM, or
#' GDINA).  \code{p0} is the lower-bound (guessing) probability and
#' \code{p1} is the upper-bound (1 - slipping) probability.
#'
#' @param patterns A binary matrix of attribute mastery patterns.
#' @param model One of \code{"DINA"}, \code{"DINO"}, \code{"ACDM"}, \code{"GDINA"}.
#' @param p0 Lower-bound probability (guessing).
#' @param p1 Upper-bound probability (1 - slipping).
#' @param mono.constraint Logical; if \code{TRUE} (default), enforce the
#'   monotonicity constraint that mastering more attributes never decreases
#'   the success probability.  Only applies to \code{model = "GDINA"}.
#' @return A named numeric vector of probabilities, one per pattern row.
#' @keywords internal
make_item_prob_cdm <- function(patterns, model, p0, p1,
                                mono.constraint = TRUE) {
  model <- toupper(model)
  patterns <- as.matrix(patterns)

  K <- ncol(patterns)
  score <- rowSums(patterns)

  if (p0 > p1) {
    stop("'P0' should be less than or equal to 'P1'.", call. = FALSE)
  }

  p <- switch(
    model,

    DINA = ifelse(score == K, p1, p0),

    DINO = ifelse(score > 0, p1, p0),

    ACDM = {
      if (K > 0) {
        w <- runif(K)
        w <- w / sum(w)
        as.vector(p0 + patterns %*% ((p1 - p0) * w))
      } else {
        rep(p0, length(score))
      }
    },

    GDINA = {
      if (K <= 1L) {
        out <- c(p0, p1)
      } else if (mono.constraint) {
        # Enforce monotonicity: mastering more attributes never
        # decreases the success probability (matching GDINA package).
        out <- rep(p0, length(score))
        out[score == K] <- p1
        mid <- which(score > 0 & score < K)
        pre <- get_precedent_patterns(patterns)
        for (l in mid) {
          out[l] <- runif(1, max(out[pre[[l]]]), p1)
        }
      } else {
        # No monotonicity: random probabilities in [p0, p1]
        out <- rep(p0, length(score))
        out[score == K] <- p1
        mid <- score > 0 & score < K
        out[mid] <- runif(sum(mid), p0, p1)
      }
      out
    }
  )

  p <- clip_prob(p)
  names(p) <- pattern_key(patterns)
  p
}


# ---- Forced-choice transformation (R implementation) --------------------------

#' Compute forced-choice pattern probabilities from item endorsement probabilities
#'
#' Transforms item-level endorsement probabilities into block-level ranking
#' pattern probabilities via the sequential Luce/Plackett choice rule.
#' This is an R implementation of the C++ \code{forced_choice_from_agree}
#' function used by \code{cpp_model_FCMIRT} and \code{cpp_model_FCGGUM}.
#'
#' @param agree An \eqn{N \times I} matrix of item endorsement probabilities.
#' @param patterns.total List of full-ranking pattern matrices (one per block).
#' @param patterns List of observed-pattern matrices (one per block).
#' @return An \eqn{N \times \sum_B P_b} matrix of forced-choice pattern
#'   probabilities, where \eqn{P_b} is the number of observed patterns for
#'   block \eqn{b}.
#' @keywords internal
forced_choice_from_agree <- function(agree, patterns.total, patterns) {
  agree <- as.matrix(agree)
  N <- nrow(agree)
  N.block <- length(patterns.total)

  if (length(patterns) != N.block) {
    stop("'patterns.total' and 'patterns' must have the same length.", call. = FALSE)
  }

  # Logit transform with numerical clamping (preserve matrix dimensions)
  eps <- 1e-12
  agree.clamped <- agree
  agree.clamped[agree.clamped < eps] <- eps
  agree.clamped[agree.clamped > 1 - eps] <- 1 - eps
  logit <- log(agree.clamped) - log1p(-agree.clamped)

  # Total number of observed pattern columns
  prob.cols <- sum(vapply(patterns, nrow, integer(1L)))
  prob <- matrix(NA_real_, N, prob.cols)

  idx <- 0L

  for (b in seq_len(N.block)) {
    pat.total <- patterns.total[[b]]
    pat.obs <- patterns[[b]]

    bR.total <- nrow(pat.total)
    bC.total <- ncol(pat.total)
    bR <- nrow(pat.obs)
    bC <- ncol(pat.obs)

    if (bR.total <= 0L || bC.total <= 0L || bR <= 0L || bC <= 0L) {
      stop("Pattern matrices must have positive dimensions.", call. = FALSE)
    }

    # Compute log-probability of each full ranking for each person
    log.prob.total <- matrix(0, N, bR.total)

    for (br in seq_len(bR.total)) {
      for (p in seq_len(N)) {
        log.prob.br <- 0.0

        for (bc in seq_len(bC.total - 1L)) {
          selected <- as.integer(pat.total[br, bc])

          # Log-sum-exp over remaining candidates
          candidates <- as.integer(pat.total[br, bc:bC.total])
          terms <- logit[p, candidates]
          max.term <- max(terms)
          denom.sum <- sum(exp(terms - max.term))
          log.denom <- max.term + log(denom.sum)

          log.prob.br <- log.prob.br + logit[p, selected] - log.denom
        }

        log.prob.total[p, br] <- log.prob.br
      }
    }

    # Softmax normalization of full ranking probabilities
    prob.total <- matrix(0, N, bR.total)
    for (p in seq_len(N)) {
      max.log <- max(log.prob.total[p, ])
      denom.sum <- sum(exp(log.prob.total[p, ] - max.log))
      log.denom <- max.log + log(denom.sum)
      prob.total[p, ] <- exp(log.prob.total[p, ] - log.denom)
    }

    # Aggregate full rankings to observed patterns
    cols <- idx + seq_len(bR)

    if (bC == bC.total) {
      # RANK: direct match; find each observed pattern in full rankings
      matches <- integer(bR)
      for (br in seq_len(bR)) {
        for (fr in seq_len(bR.total)) {
          if (all(pat.obs[br, ] == pat.total[fr, ])) {
            matches[br] <- fr
            break
          }
        }
        if (matches[br] <= 0L) {
          stop("'patterns' contains a full ranking absent from 'patterns.total'.", call. = FALSE)
        }
      }
      for (br in seq_len(bR)) {
        prob[, cols[br]] <- prob.total[, matches[br]]
      }

    } else if (bC == 2L) {
      # MOLE: sum full rankings where first matches AND last matches
      matches <- vector("list", bR)
      for (br in seq_len(bR)) {
        first <- as.integer(pat.obs[br, 1L])
        last  <- as.integer(pat.obs[br, 2L])
        matched.fr <- integer(0L)
        for (fr in seq_len(bR.total)) {
          full.first <- as.integer(pat.total[fr, 1L])
          full.last  <- as.integer(pat.total[fr, bC.total])
          if (full.first == first && full.last == last) {
            matched.fr <- c(matched.fr, fr)
          }
        }
        matches[[br]] <- matched.fr
      }
      for (br in seq_len(bR)) {
        if (length(matches[[br]]) > 0L) {
          prob[, cols[br]] <- rowSums(
            prob.total[, matches[[br]], drop = FALSE]
          )
        } else {
          prob[, cols[br]] <- 0.0
        }
      }

    } else if (bC == 1L) {
      # PICK: sum full rankings where first item matches
      matches <- vector("list", bR)
      for (br in seq_len(bR)) {
        first <- as.integer(pat.obs[br, 1L])
        matched.fr <- integer(0L)
        for (fr in seq_len(bR.total)) {
          full.first <- as.integer(pat.total[fr, 1L])
          if (full.first == first) {
            matched.fr <- c(matched.fr, fr)
          }
        }
        matches[[br]] <- matched.fr
      }
      for (br in seq_len(bR)) {
        if (length(matches[[br]]) > 0L) {
          prob[, cols[br]] <- rowSums(
            prob.total[, matches[[br]], drop = FALSE]
          )
        } else {
          prob[, cols[br]] <- 0.0
        }
      }

    } else {
      stop("Unsupported forced-choice pattern length.", call. = FALSE)
    }

    # Normalize observed pattern probabilities within each block
    for (p in seq_len(N)) {
      row.sum <- sum(prob[p, cols])
      if (row.sum <= 0.0 || !is.finite(row.sum)) {
        stop("Computed forced-choice probabilities cannot be normalized.", call. = FALSE)
      }
      prob[p, cols] <- prob[p, cols] / row.sum
    }

    idx <- idx + bR
  }

  return(prob)
}


# ---- CDM monotonicity constraint ----------------------------------------------

#' Find precedent (subset) patterns for monotonicity enforcement
#'
#' For each attribute mastery pattern, identifies which other patterns are
#' proper subsets (i.e., mastering a subset of the attributes).  Used to
#' enforce the GDINA monotonicity constraint: mastering more attributes
#' should never decrease the success probability.
#'
#' Equivalent to \code{preloclist} in the GDINA R package.
#'
#' @param patterns A binary matrix of attribute mastery patterns
#'   (\eqn{L \times K}).
#' @return A list of length \eqn{L}; the \eqn{l}-th element is an integer
#'   vector of row indices corresponding to patterns that are proper
#'   subsets of pattern \eqn{l}.
#' @keywords internal
get_precedent_patterns <- function(patterns) {
  patterns <- as.matrix(patterns)
  L <- nrow(patterns)
  pre <- vector("list", L)
  for (l in seq_len(L)) {
    x <- patterns[l, ]
    # A pattern j is a subset of pattern l if:
    #   x OR patterns[j,] == x   AND   patterns[j,] != x
    is_subset <- rep(FALSE, L)
    for (j in seq_len(L)) {
      if (j == l) next
      if (all((x | patterns[j, ]) == x) && any(patterns[j, ] != x)) {
        is_subset[j] <- TRUE
      }
    }
    pre[[l]] <- which(is_subset)
  }
  pre
}

# ---- CDM design matrix --------------------------------------------------------

#' Build the CDM design matrix for a set of attribute patterns
#'
#' Constructs the design matrix \eqn{X} used in the cognitive diagnostic
#' model probability formulation \eqn{p = X\delta}.  The structure depends
#' on the model type:
#' \describe{
#'   \item{DINA}{Intercept + single mastery indicator (all attributes mastered).}
#'   \item{DINO}{Intercept + single "at-least-one" indicator.}
#'   \item{ACDM}{Intercept + main effects of each attribute.}
#'   \item{GDINA}{Full saturated interaction matrix (all \eqn{2^K} columns).}
#' }
#'
#' @param patterns A binary matrix of attribute mastery patterns
#'   (\eqn{L \times K}).
#' @param model One of \code{"DINA"}, \code{"DINO"}, \code{"ACDM"}, \code{"GDINA"}.
#' @return A numeric design matrix \eqn{X} with \eqn{L} rows.  Column count
#'   depends on the model.
#' @keywords internal
get_design_matrix_cdm <- function(patterns, model = "GDINA") {
  model <- toupper(model)
  patterns <- as.matrix(patterns)

  L <- nrow(patterns)
  K <- ncol(patterns)

  if (!model %in% c("DINA", "DINO", "ACDM", "GDINA")) {
    stop("'model' must be one of: DINA, DINO, ACDM, GDINA.", call. = FALSE)
  }

  if (model == "DINA") {
    X <- cbind(
      intercept = 1,
      master    = as.integer(rowSums(patterns) == K)
    )
  } else if (model == "DINO") {
    X <- cbind(
      intercept    = 1,
      at_least_one = as.integer(rowSums(patterns) > 0)
    )
  } else if (model == "ACDM") {
    X <- cbind(intercept = 1, patterns)
    colnames(X) <- c("intercept", paste0("main", seq_len(K)))
  } else {
    # GDINA: full interaction matrix
    X <- outer(
      seq_len(L),
      seq_len(L),
      Vectorize(function(l1, l2) prod(patterns[l1, ] ^ patterns[l2, ]))
    )
    colnames(X) <- paste0("d", pattern_key(patterns))
  }

  storage.mode(X) <- "numeric"
  X
}
