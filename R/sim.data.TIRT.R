#' Simulate Data from the Thurstonian IRT (TIRT) Model
#'
#' @description
#' Generates forced-choice pairwise comparison data under the Thurstonian
#' IRT framework. Each statement loads on exactly one latent dimension;
#' pairwise latent differences with probit link determine binary outcomes.
#' Supports RANK, MOLE, and PICK response types. For each pair \eqn{(i,k)},
#' the simulated comparison uses
#' \eqn{-\gamma_{ik} + \lambda_i q_i'\theta_j -
#' \lambda_k q_k'\theta_j + \epsilon_i - \epsilon_k}, with
#' \eqn{\epsilon_i \sim N(0,\psi_i^2)}.
#'
#' @param N.person Integer; number of persons (default: 1000).
#' @param N.block Integer; number of blocks (default: 10).
#' @param I.block Integer; items per block (default: 2).
#' @param D Integer; number of latent dimensions (default: 3, >= 2).
#' @param fc.type Character; \code{"RANK"} (default), \code{"MOLE"}, or
#'   \code{"PICK"}.
#' @param control Optional list with entries \code{Q.matrix}, \code{block.items},
#'   \code{lambda}, \code{psi2}, \code{gamma.matrix}, \code{gamma.item},
#'   \code{Corr}.
#'
#' @return An object of class \code{"data.TIRT"} containing \code{data},
#'   \code{response}, \code{theta}, \code{par}, \code{Q.matrix},
#'   \code{block.items}, \code{Corr}, \code{gamma.matrix}, \code{lambda},
#'   \code{psi2}, \code{pairs.value}, and data-generating arguments.
#'
#' @seealso \code{\link{fit.TIRT}}, \code{\link{model.TIRT}}
#'
#' @examples
#' set.seed(123)
#' sim <- sim.data.TIRT(N.person = 20, N.block = 3, I.block = 2,
#'                      D = 2, fc.type = "RANK")
#' str(sim$data)
#' head(sim$response)
#' dim(sim$Q.matrix)
#'
#' @export
sim.data.TIRT <- function(N.person = 1000, N.block = 10, I.block = 2, D = 3,
                          fc.type = "RANK", control = NULL) {

  call <- match.call()

  if (D < 2) {
    stop("TIRT requires at least 2 latent dimensions (D >= 2).", call. = FALSE)
  }

  if (is.null(control)) control <- list()

  # ---- Scalar integer validation ----
  N.person <- check_integer_scalar(N.person, "N.person", 1L)
  D        <- check_integer_scalar(D,        "D",        2L)
  I.block  <- check_integer_scalar(I.block,  "I.block",  2L)
  N.block  <- check_integer_scalar(N.block,  "N.block",  1L)

  # ---- Q-matrix and block configuration ----
  Q.matrix    <- control$Q.matrix
  block.items <- control$block.items

  res <- resolve_fc_blocks(
    Q.matrix       = Q.matrix,
    block.items    = block.items,
    I.block        = I.block,
    N.block        = N.block,
    D              = D,
    q.valid.values = c(-1, 0, 1),
    q.row.check    = function(r) sum(abs(r)) == 1L
  )
  I.states   <- res$I.states
  N.block    <- res$N.block
  I.block    <- res$I.block
  block.items <- res$block.items
  Q.matrix   <- res$Q.matrix
  D          <- res$D

  if (D < 2) {
    stop("TIRT requires at least 2 latent dimensions (D >= 2).", call. = FALSE)
  }

  # ---- fc.type validation ----
  fc.type <- normalize_fc_type(fc.type, N.block)

  # ---- Corr and theta ----
  Corr  <- validate_corr_matrix(control$Corr, D)
  theta <- generate_theta_mvn(N.person, D, Corr)

  # ---- Auto-generate Q-matrix if not supplied ----
  if (is.null(Q.matrix)) {
    Q.matrix <- sim.Q.matrix.FC(I.states = I.states, D = D,
                                block.items = block.items,
                                allow.negative = TRUE)
  }

  # ---- TIRT-specific: lambda, psi2, gamma ----
  lambda <- control$lambda
  psi2   <- control$psi2
  gamma.matrix <- control$gamma.matrix

  last_item_in_block <- as.integer(vapply(block.items, function(x) x[length(x)], numeric(1L)))

  if (is.null(lambda)) {
    lambda <- runif(I.states, 0.40, 0.90)
    if (D == 2 && length(I.block) == 1L && I.block == 2) {
      first_pair_items <- block.items[[1]]
      lambda[first_pair_items] <- rep(0.8, length(first_pair_items))
    }
  } else {
    if (length(lambda) != I.states) stop("Length of 'lambda' must be ", I.states, call. = FALSE)
    if (D == 2 && length(I.block) == 1L && I.block == 2) {
      first_pair_items <- block.items[[1]]
      lambda[first_pair_items] <- rep(0.8, length(first_pair_items))
    }
  }

  is_case_a <- (length(I.block) == 1L && I.block == 2 && D > 2)
  is_case_b <- (D == 2 && length(I.block) == 1L && I.block == 2)

  if (is.null(psi2)) {
    psi2 <- 1 - lambda^2
    if (!is_case_a && !is_case_b) {
      psi2[last_item_in_block] <- 1.0
    }
    if (is_case_a || is_case_b) {
      psi2[] <- 0.5
    }
  } else {
    if (length(psi2) != I.states) stop("Length of 'psi2' must be ", I.states, call. = FALSE)
    if (!is_case_a && !is_case_b) {
      psi2[last_item_in_block] <- 1.0
    }
    if (is_case_a || is_case_b) {
      psi2[] <- 0.5
    }
  }

  if (is.null(gamma.matrix)) {
    gamma.matrix <- matrix(NA, I.states, I.states)
  } else {
    if (!is.matrix(gamma.matrix) || nrow(gamma.matrix) != I.states ||
        ncol(gamma.matrix) != I.states)
      stop("'gamma.matrix' must be an ", I.states, "x", I.states, " matrix", call. = FALSE)
  }

  # ---- TIRT-specific: gamma generation ----
  # Generate item-level gamma utilities to ensure transitive pairwise comparisons
  if (all(is.na(gamma.matrix))) {
    gamma.item <- get_ctrl("gamma.item", runif(I.states, -1, 1), control)
    if (length(gamma.item) != I.states)
      stop("'gamma.item' must have length ", I.states, call. = FALSE)
    for (i in seq_len(I.states)) {
      for (k in seq_len(I.states)) {
        if (i != k) {
          gamma.matrix[i, k] <- gamma.item[i] - gamma.item[k]
        }
      }
    }
  } else {
    # User-provided gamma.matrix: check for intransitivities and warn
    for (i in seq_len(I.states)) {
      for (j in seq_len(I.states)) {
        for (k in seq_len(I.states)) {
          if (i != j && j != k && i != k) {
            if (!is.na(gamma.matrix[i, j]) && !is.na(gamma.matrix[j, k]) &&
                !is.na(gamma.matrix[i, k])) {
              sum_ij_jk <- gamma.matrix[i, j] + gamma.matrix[j, k]
              if (abs(sum_ij_jk - gamma.matrix[i, k]) > 1e-6) {
                warning(sprintf(
                  "Intransitivity detected in gamma.matrix: gamma[%d,%d]=%.4f + gamma[%d,%d]=%.4f = %.4f != gamma[%d,%d]=%.4f.",
                  i, j, gamma.matrix[i, j], j, k, gamma.matrix[j, k],
                  sum_ij_jk, i, k, gamma.matrix[i, k]
                ))
              }
            }
          }
        }
      }
    }
  }

  # ---- TIRT-specific: pairwise response simulation ----
  pairs.total <- sum(vapply(block.items, function(x) choose(length(x), 2), numeric(1L)))
  response.total <- matrix(NA, N.person, pairs.total)
  pairs <- c()

  gamma <- rep(NA, pairs.total)
  idx <- 1L
  for (b in seq_len(N.block)) {
    block.items.cur <- block.items[[b]]
    pairs.value.total <- t(combn(block.items.cur, 2))
    I.pairs.cur <- nrow(pairs.value.total)
    for (pp in seq_len(I.pairs.cur)) {
      i <- pairs.value.total[pp, 1]
      k <- pairs.value.total[pp, 2]
      if (is.na(gamma.matrix[i, k])) {
        if (!exists("gamma.item", inherits = FALSE)) {
          gamma.item <- runif(I.states, -1.0, 1.0)
        }
        gamma_val <- gamma.item[i] - gamma.item[k]
        gamma.matrix[i, k] <- gamma_val
        gamma.matrix[k, i] <- -gamma_val
      }
      gamma[idx] <- gamma.matrix[i, k]
      idx <- idx + 1L
    }

    if (fc.type[b] == "RANK") {
      pairs <- c(pairs, I.pairs.cur)
    } else if (fc.type[b] == "MOLE") {
      pairs <- c(pairs, length(block.items.cur) * 2 - 3)
    } else if (fc.type[b] == "PICK") {
      pairs <- c(pairs, length(block.items.cur) - 1)
    }
  }

  response <- matrix(NA, N.person, sum(pairs))
  pairs.value <- vector("list", N.person)
  for (p in seq_len(N.person)) {
    item_errors <- rnorm(I.states, mean = 0, sd = sqrt(psi2))
    idx.total <- 0L
    idx <- 0L
    for (b in seq_len(N.block)) {
      block.items.cur <- block.items[[b]]
      pairs.value.total <- t(combn(block.items.cur, 2))

      I.pairs.total.cur <- nrow(pairs.value.total)
      scores.b <- rep(0, length(block.items.cur))
      for (pp in seq_len(I.pairs.total.cur)) {
        i <- pairs.value.total[pp, 1]
        k <- pairs.value.total[pp, 2]

        lambda.i <- lambda[i]
        lambda.k <- lambda[k]
        gamma.current <- gamma.matrix[i, k]

        mu <- -gamma.current +
          lambda.i * sum(Q.matrix[i, ] * theta[p, ]) -
          lambda.k * sum(Q.matrix[k, ] * theta[p, ])
        latent_y <- mu + item_errors[i] - item_errors[k]
        response.total[p, idx.total + pp] <- as.integer(latent_y >= 0)

        if (response.total[p, idx.total + pp] > 0) {
          scores.b[match(i, block.items.cur)] <- scores.b[match(i, block.items.cur)] + 1L
        } else {
          scores.b[match(k, block.items.cur)] <- scores.b[match(k, block.items.cur)] + 1L
        }
      }

      if (fc.type[b] == "RANK") {
        pairs.value[[p]][[b]] <- pairs.value.total
        response[p, (idx + 1):(idx + pairs[b])] <-
          response.total[p, (idx.total + 1):(idx.total + I.pairs.total.cur)]
      } else if (fc.type[b] == "MOLE") {
        pairs.max <- which.max(scores.b)[1]
        pairs.min <- which.min(scores.b)[1]

        pairs.value.posi <- sort(unique(c(
          seq_len(I.pairs.total.cur) * (pairs.value.total[, 1] == block.items.cur[pairs.max]),
          seq_len(I.pairs.total.cur) * (pairs.value.total[, 1] == block.items.cur[pairs.min]),
          seq_len(I.pairs.total.cur) * (pairs.value.total[, 2] == block.items.cur[pairs.max]),
          seq_len(I.pairs.total.cur) * (pairs.value.total[, 2] == block.items.cur[pairs.min])
        )))
        pairs.value.posi <- pairs.value.posi[pairs.value.posi != 0]
        pairs.value[[p]][[b]] <- pairs.value.total[pairs.value.posi, , drop = FALSE]

        response[p, (idx + 1):(idx + pairs[b])] <-
          response.total[p, (idx.total + seq_len(I.pairs.total.cur))[pairs.value.posi]]
      } else if (fc.type[b] == "PICK") {
        pairs.max <- which.max(scores.b)[1]
        pairs.value.posi <- sort(unique(c(
          seq_len(I.pairs.total.cur) * (pairs.value.total[, 1] == block.items.cur[pairs.max]),
          seq_len(I.pairs.total.cur) * (pairs.value.total[, 2] == block.items.cur[pairs.max])
        )))
        pairs.value.posi <- pairs.value.posi[pairs.value.posi != 0]
        pairs.value[[p]][[b]] <- pairs.value.total[pairs.value.posi, , drop = FALSE]

        response[p, (idx + 1):(idx + pairs[b])] <-
          response.total[p, (idx.total + seq_len(I.pairs.total.cur))[pairs.value.posi]]
      }
      idx.total <- idx.total + I.pairs.total.cur
      idx <- idx + pairs[b]
    }
  }

  # ---- Row/column names ----
  colnames(response) <- paste0("pairs", seq_len(sum(pairs)))
  rownames(Q.matrix) <- rownames(gamma.matrix) <-
    colnames(gamma.matrix) <- paste0("S", seq_len(I.states))
  colnames(Q.matrix) <- colnames(theta) <- paste0("Dim.", seq_len(D))

  lambda_sign <- apply(Q.matrix, 1, function(x) {
    nz <- which(x != 0)
    if (length(nz) == 0) return(1)
    sign(x[nz[1]])
  })
  lambda_signed <- lambda * lambda_sign
  par <- cbind(lambda_signed * abs(Q.matrix), psi2)

  data <- get.data.from.response.TIRT(response, block.items, fc.type, pairs.value)

  # ---- Output assembly ----
  results <- list(
    data         = data,
    response     = response,
    theta        = theta,
    par          = par,
    Q.matrix     = Q.matrix,
    block.items  = block.items,
    Corr         = Corr,
    pairs.value  = pairs.value,
    gamma.matrix = gamma.matrix,
    gamma        = gamma,
    gamma.item   = if (exists("gamma.item", inherits = FALSE)) gamma.item else NULL,
    lambda       = lambda_signed,
    psi2         = psi2,
    N.person     = N.person,
    N.block      = N.block,
    I.block      = I.block,
    D            = D,
    I.states     = I.states,
    fc.type      = fc.type,
    call         = call,
    arguments    = list(
      N.person = N.person,
      N.block  = N.block,
      I.block  = I.block,
      D        = D,
      fc.type  = fc.type,
      control  = control
    )
  )

  class(results) <- "data.TIRT"
  return(results)
}

#' Compute TIRT Pairwise Comparison Probabilities
#' @param theta Latent trait matrix.
#' @param lambda Factor loading magnitudes per statement. If signed values
#'   are supplied from a fitted object, the function uses \code{abs(lambda)}
#'   and the sign in \code{Q.matrix} carries the direction.
#' @param psi2 Uniqueness variances per statement.
#' @param gamma.matrix Pairwise intercept matrix.
#' @param Q.matrix Design Q-matrix.
#' @param pairs.matrix Matrix of statement-index pairs.
#' @return Matrix of pairwise choice probabilities:
#'   \eqn{\Phi\{[-\gamma_{ik} + (|\lambda_i|q_i -
#'   |\lambda_k|q_k)'\theta_j] / \sqrt{\psi_i^2+\psi_k^2}\}}.
#' @keywords internal
model.TIRT <- function(theta, lambda, psi2, gamma.matrix, Q.matrix, pairs.matrix) {
  lambda.matrix <- t(abs(lambda) * Q.matrix)
  N.person <- nrow(theta)
  I.pairs <- nrow(pairs.matrix)
  prob <- matrix(0, N.person, I.pairs)

  for (i in seq_len(I.pairs)) {
    mu <- -gamma.matrix[pairs.matrix[i, 1], pairs.matrix[i, 2]] +
      theta %*% lambda.matrix[, pairs.matrix[i, 1]] -
      theta %*% lambda.matrix[, pairs.matrix[i, 2]]
    pis2l <- sqrt(max(psi2[pairs.matrix[i, 1]] + psi2[pairs.matrix[i, 2]],
                      .Machine$double.eps))
    prob[, i] <- pnorm(mu / pis2l)
  }
  return(prob)
}
