#' Simulate Data from the Forced-Choice GGUM Model
#'
#' @description
#' Generates forced-choice ranking data using the GGUM ideal-point model as the
#' item-level endorsement process. Within-block rankings follow the same
#' sequential Luce/Plackett mechanism as \code{\link{model.FCGGUM}}: the
#' utility for each statement is the logit of its binary GGUM endorsement
#' probability.
#'
#' @param N.person Integer; number of persons (default: 1000).
#' @param N.block Integer; number of forced-choice blocks (default: 10).
#' @param I.block Integer; items per block (default: 2).
#' @param D Integer; number of latent dimensions (default: 3, must be >= 2).
#' @param fc.type Character; \code{"RANK"} (default), \code{"MOLE"}, or
#'   \code{"PICK"}.
#' @param control Optional list with entries \code{Q.matrix},
#'   \code{block.items}, and \code{Corr}.
#'
#' @return An object of class \code{"data.FCGGUM"} containing \code{data},
#'   \code{response}, \code{theta}, \code{par}, \code{Q.matrix},
#'   \code{block.items}, \code{Corr}, \code{patterns}, \code{patterns.total},
#'   \code{prob}, \code{mask}, and data-generating arguments.
#'
#' @seealso \code{\link{fit.FCGGUM}}, \code{\link{model.FCGGUM}}
#'
#' @examples
#' set.seed(123)
#' sim <- sim.data.FCGGUM(N.person = 20, N.block = 3, I.block = 2,
#'                        D = 2, fc.type = "RANK")
#' str(sim$data)
#' head(sim$response)
#' dim(sim$prob)
#'
#' @export
sim.data.FCGGUM <- function(N.person = 1000, N.block = 10, I.block = 2,
                            D = 3, fc.type = "RANK", control = NULL) {

  call <- match.call()

  if (D < 2) {
    stop("FCGGUM requires at least 2 latent dimensions (D >= 2).", call. = FALSE)
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
    q.row.check    = function(r) any(r != 0)
  )
  I.states   <- res$I.states
  N.block    <- res$N.block
  I.block    <- res$I.block
  block.items <- res$block.items
  Q.matrix   <- res$Q.matrix
  D          <- res$D

  if (D < 2) {
    stop("FCGGUM requires at least 2 latent dimensions (D >= 2).", call. = FALSE)
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

  # ---- GGUM item parameters ----
  a <- matrix(runif(I.states * D, 0.5, 2.0), I.states, D) * abs(Q.matrix)

  par <- matrix(NA, nrow = I.states, ncol = D + D + 2)
  par[, 1:D] <- a
  mask <- matrix(NA, nrow = I.states, ncol = D + D + 2)
  mask[, 1:D] <- 1
  for (i in seq_len(I.states)) {

    for (d in seq_len(D)) {
      if (Q.matrix[i, d] == 1) {
        par[i, D + d] <- runif(1, 0.0, 2.0)
      } else if (Q.matrix[i, d] == -1) {
        par[i, D + d] <- runif(1, -2.0, 0.0)
      } else {
        par[i, D + d] <- 0.0
      }
    }
    mask[i, (D + 1):(D + D)] <- 1

    par[i, D + D + 1] <- 0.0
    mask[i, D + D + 1] <- 1
    Ki <- 2
    Kf <- Ki - 1
    tau_free <- numeric(Kf)
    for (t in seq_len(Kf)) {
      lo <- -2.00 + (t - 1) * 2.0 / Kf
      hi <- -2.00 +  t      * 2.0 / Kf
      tau_free[t] <- runif(1, lo, hi)
    }

    tau_sorted <- sort(tau_free, decreasing = FALSE)
    par[i, (D + D + 2):(D + D + Ki)] <- tau_sorted
    mask[i, (D + D + 2):(D + D + Ki)] <- 1
  }

  # ---- FC permutation patterns ----
  pat <- generate_fc_permutation_patterns(block.items, fc.type)
  patterns.total <- pat$patterns.total
  patterns       <- pat$patterns

  # ---- FC pattern probabilities ----
  prob <- model.FCGGUM(theta, par, patterns.total, patterns)

  # ---- Response sampling ----
  samples <- sample_fc_responses(prob, patterns, N.person, N.block)
  response <- samples$response
  data     <- samples$data

  # ---- Row/column names ----
  rownames(theta) <- paste0("person", seq_len(N.person))
  if (D > 1) {
    colnames(theta) <- paste0("theta", seq_len(D))
  }

  rownames(par) <- paste0("item", seq_len(I.states))
  colnames(par) <- c(paste0("a", seq_len(D)),
                     paste0("delta", seq_len(D)),
                     paste0("tau", 0:1))
  rownames(Q.matrix) <- paste0("item", seq_len(I.states))
  colnames(Q.matrix) <- paste0("dim", seq_len(D))

  # ---- Output assembly ----
  data.obj <- list(
    data           = data,
    response       = response,
    theta          = theta,
    par            = par,
    Q.matrix       = Q.matrix,
    block.items    = block.items,
    Corr           = Corr,
    patterns       = patterns,
    patterns.total = patterns.total,
    prob           = prob,
    mask           = mask,
    N.person       = N.person,
    N.block        = N.block,
    I.block        = I.block,
    D              = D,
    I.states       = I.states,
    fc.type        = fc.type,
    call           = call,
    arguments      = list(
      N.person = N.person,
      N.block  = N.block,
      I.block  = I.block,
      D        = D,
      fc.type  = fc.type,
      control  = control
    )
  )

  class(data.obj) <- "data.FCGGUM"

  return(data.obj)
}

#' Compute FCGGUM Forced-Choice Pattern Probabilities
#' @param theta Latent trait matrix.
#' @param par Binary GGUM item parameter matrix with columns
#'   \code{a1..aD, delta1..deltaD, tau0, tau1}; \code{tau0} is fixed at 0
#'   and \code{tau1 < 0}.
#' @param patterns.total Full ranking patterns.
#' @param patterns Observed ranking patterns.
#' @return Block ranking probability matrix. The code first computes binary
#'   GGUM endorsement probabilities from the distance-based MGGUM formula,
#'   transforms them to logits, and then applies the sequential
#'   Luce/Plackett rule. MOLE and PICK probabilities are sums over
#'   compatible full rankings normalized over observed patterns.
#' @keywords internal
model.FCGGUM <- function(theta, par, patterns.total, patterns){
  cpp_model_FCGGUM(theta, par, patterns.total, patterns)
}
