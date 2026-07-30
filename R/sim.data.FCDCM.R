#' Simulate Data from the Forced-Choice DCM Model
#'
#' @description
#' Generates forced-choice paired-comparison data from a higher-order
#' cognitive diagnostic model (DINA/DINO). The continuous higher-order
#' trait governs attribute mastery probabilities via a logistic link;
#' two-statement blocks produce binary choices based on attribute
#' mastery patterns. In each block, selecting the first statement has
#' probability \eqn{\eta_0} when its condensed mastery is lower than the
#' second statement's, \eqn{0.5} when they are equal, and
#' \eqn{0.5+\eta_{AB}} when it is higher.
#'
#' @param N.person Integer; number of persons (default: 1000).
#' @param N.block Integer; number of forced-choice block pairs (default: 10).
#' @param D Integer; number of attributes (default: 3, must be >= 2).
#' @param dcm.type Character; \code{"DINA"} (default) or \code{"DINO"}.
#' @param control Optional list with entries \code{Q.matrix}, \code{block.items},
#'   \code{single}, \code{delta1}, \code{delta0},
#'   \code{theta}, \code{alpha}, \code{par}.  If \code{block.items} is
#'   omitted, statements are assembled into cross-Q-vector forced-choice
#'   blocks with balanced first/second positions and diversified attribute
#'   pairs.  \code{par} is a \eqn{B \times 2} matrix with columns
#'   \code{eta0} and \code{etaAB}.
#'
#' @return An object of class \code{"data.FCDCM"} containing \code{data},
#'   \code{response}, \code{theta}, \code{alpha}, \code{zeta}, \code{par},
#'   \code{Q.matrix}, \code{block.items},
#'   \code{prob}, and data-generating arguments.
#'
#' @seealso \code{\link{fit.FCDCM}}, \code{\link{model.FCDCM}}
#'
#' @examples
#' set.seed(123)
#' sim <- sim.data.FCDCM(N.person = 20, N.block = 3, D = 2,
#'                       dcm.type = "DINA")
#' str(sim$data)
#' head(sim$response)
#' head(sim$alpha)
#'
#' @export
sim.data.FCDCM <- function(N.person = 1000, N.block = 10, D = 3,
                           dcm.type = "DINA", control = NULL) {

  call <- match.call()

  if (is.null(control)) {
    control <- list()
  }

  N.person <- as.integer(N.person[1L])
  N.block <- as.integer(N.block[1L])
  D <- as.integer(D[1L])
  if (is.na(N.person) || N.person < 1L ||
      is.na(N.block) || N.block < 1L ||
      is.na(D) || D < 2L) {
    stop("'N.person', 'N.block', and 'D' must be valid positive integers.", call. = FALSE)
  }
  I.states <- 2 * N.block

  Q.matrix    <- control$Q.matrix
  single      <- control$single
  if (is.null(Q.matrix)) {
    if (is.null(single)) {
      single <- TRUE
    }
    Q.matrix <- sim.data.Q.CDM(D, I.states, single = single)
  }

  block.items <- control$block.items
  if(is.null(block.items)){
    # Group items by their Q-vector pattern and interleave to ensure
    # each block pairs items with DIFFERENT Q-vectors.  Blocks whose
    # two items share the same Q-vector are completely uninformative
    # under DINA/DINO (za = zb for every attribute class).
    q.hash <- apply(Q.matrix, 1L, paste, collapse = ":")
    q.groups <- split(seq_len(I.states), q.hash)
    q.groups <- lapply(q.groups, function(x) sample(x, length(x)))
    n.groups <- length(q.groups)
    if (n.groups >= 2L) {
      block.items <- fcdcm_balanced_q_blocks(q.groups, N.block)
    } else {
      interleaved <- sample(seq_len(I.states), I.states)
      block.items <- split(interleaved, rep(seq_len(N.block), each = 2L))
    }
  }
  patterns <- matrix(unlist(block.items), ncol=2, byrow=TRUE)

  dcm.type <- toupper(dcm.type)
  if (length(dcm.type) != I.states) {
    dcm.type <- rep(dcm.type[1L], I.states)
  }
  if (!all(dcm.type %in% c("DINA", "DINO"))) {
    stop("'dcm.type' must contain only 'DINA' or 'DINO'.", call. = FALSE)
  }

  delta1 <- control$delta1
  delta0 <- control$delta0
  theta  <- control$theta
  alpha  <- control$alpha

  if (is.null(alpha)) {
    if (is.null(delta1)) {
      delta1 <- rlnorm(D, meanlog = 0.25, sdlog = 0.25)
    } else {
      delta1 <- as.numeric(delta1)
      if (length(delta1) != D) {
        stop("'control$delta1' must have length ", D, ".", call. = FALSE)
      }
    }

    if (is.null(delta0)) {
      delta0 <- rnorm(D, mean = 0, sd = 1)
    } else {
      delta0 <- as.numeric(delta0)
      if (length(delta0) != D) {
        stop("'control$delta0' must have length ", D, ".", call. = FALSE)
      }
    }

    if (is.null(theta)) {
      theta <- rnorm(N.person, mean = 0, sd = 1)
    } else {
      theta <- as.numeric(theta)
      if (length(theta) != N.person) {
        stop("'control$theta' must have length ", N.person, ".", call. = FALSE)
      }
    }

    prob.alpha <- distribution_higher_order(theta = theta, delta1 = delta1, delta0 = delta0)
    alpha <- (prob.alpha > matrix(runif(N.person * D), N.person, D)) * 1L
  } else {
    alpha <- as.matrix(alpha)
    if (nrow(alpha) != N.person || ncol(alpha) != D) {
      stop("'control$alpha' must be a ", N.person, "x", D, " binary matrix.", call. = FALSE)
    }
    storage.mode(alpha) <- "integer"
    if (is.null(theta)) theta <- rep(NA_real_, N.person)
    if (is.null(delta1)) delta1 <- rep(NA_real_, D)
    if (is.null(delta0)) delta0 <- rep(NA_real_, D)
    prob.alpha <- if (all(is.finite(theta)) &&
                      all(is.finite(delta1)) &&
                      all(is.finite(delta0))) {
      distribution_higher_order(theta = theta, delta1 = delta1, delta0 = delta0)
    } else {
      matrix(NA_real_, N.person, D)
    }
  }

  zeta <- get.zeta(alpha = alpha, Q.matrix = Q.matrix, dcm.type = dcm.type)

  par <- control$par
  if(is.null(par)){
    par <- cbind(
      eta0 = runif(N.block, 0.05, 0.15),
      etaAB = runif(N.block, 0.30, 0.50)
    )
  }
  par <- fcdcm_normalize_par(par, B = N.block, arg = "control$par")

  prob <- model.FCDCM(zeta = zeta, par = par, patterns = patterns)
  response <- (prob > matrix(runif(N.person * N.block), N.person, N.block)) * 1L
  data <- as.matrix(get.data.from.response.FCDCM(response, block.items))

  block.names <- paste0("B", seq_len(N.block))
  item.names <- rownames(Q.matrix)
  if (is.null(item.names)) {
    item.names <- paste0("S", seq_len(I.states))
  }
  dim.names <- colnames(Q.matrix)
  if (is.null(dim.names)) {
    dim.names <- paste0("Dim.", seq_len(D))
  }

  names(delta1) <- names(delta0) <- dim.names


  colnames(alpha) <- colnames(prob.alpha) <- colnames(Q.matrix) <- dim.names
  colnames(zeta) <- rownames(Q.matrix) <- item.names
  colnames(response) <- colnames(data) <- colnames(prob) <- rownames(patterns) <-
    rownames(par) <- names(block.items) <- block.names
  colnames(patterns) <- c("SA", "SB")

  data.obj <- list(
    data        = data,
    response    = response,
    theta       = theta,
    alpha       = alpha,
    prob.alpha  = prob.alpha,
    zeta        = zeta,
    par         = par,
    delta1      = delta1,
    delta0      = delta0,
    Q.matrix    = Q.matrix,
    block.items = block.items,
    patterns    = patterns,
    dcm.type    = dcm.type,
    prob        = prob,
    N.person    = N.person,
    N.block     = N.block,
    I.block     = 2L,
    D           = D,
    I.states    = I.states,
    call        = call,
    arguments   = list(
      N.person = N.person,
      N.block  = N.block,
      D        = D,
      dcm.type = dcm.type,
      control  = control
    )
  )

  class(data.obj) <- "data.FCDCM"

  return(data.obj)
}

fcdcm_balanced_q_blocks <- function(q.groups, N.block) {
  q.keys <- names(q.groups)
  n.item <- vapply(q.groups, length, integer(1L))
  if (length(q.keys) == 2L && n.item[1L] == n.item[2L] &&
      n.item[1L] %% 2L == 0L) {
    n.first <- n.item[1L] %/% 2L
    blocks <- vector("list", N.block)
    pos <- 1L
    for (i in seq_len(n.first)) {
      blocks[[pos]] <- c(q.groups[[q.keys[1L]]][i],
                         q.groups[[q.keys[2L]]][n.first + i])
      pos <- pos + 1L
      blocks[[pos]] <- c(q.groups[[q.keys[2L]]][i],
                         q.groups[[q.keys[1L]]][n.first + i])
      pos <- pos + 1L
    }
    return(blocks[sample(seq_along(blocks))])
  }

  first.n <- n.item %/% 2L
  extra <- N.block - sum(first.n)
  if (extra > 0L) {
    add <- order(n.item %% 2L, n.item, decreasing = TRUE)[seq_len(extra)]
    first.n[add] <- first.n[add] + 1L
  }
  second.n <- n.item - first.n

  first.pool <- second.pool <- vector("list", length(q.keys))
  names(first.pool) <- names(second.pool) <- q.keys
  for (key in q.keys) {
    first.pool[[key]] <- q.groups[[key]][seq_len(first.n[key])]
    second.pool[[key]] <- q.groups[[key]][first.n[key] + seq_len(second.n[key])]
  }

  pair.keys <- combn(q.keys, 2L, function(x) paste(sort(x), collapse = "\r"))
  ordered.grid <- expand.grid(first = q.keys, second = q.keys,
                              stringsAsFactors = FALSE)
  ordered.grid <- ordered.grid[ordered.grid$first != ordered.grid$second, ,
                               drop = FALSE]
  ordered.keys <- paste(ordered.grid$first, ordered.grid$second, sep = "\r")

  build_once <- function() {
    fp <- first.pool
    sp <- second.pool
    pair.count <- setNames(integer(length(pair.keys)), pair.keys)
    ordered.count <- setNames(integer(length(ordered.keys)), ordered.keys)
    out <- vector("list", N.block)
    key.out <- matrix(NA_character_, N.block, 2L)

    for (b in seq_len(N.block)) {
      first.left <- vapply(fp, length, integer(1L))
      second.left <- vapply(sp, length, integer(1L))
      first.available <- q.keys[first.left > 0L]
      second.available <- q.keys[second.left > 0L]
      candidates <- expand.grid(first = first.available,
                                second = second.available,
                                stringsAsFactors = FALSE)
      score <- vapply(seq_len(nrow(candidates)), function(i) {
        first.key <- candidates$first[i]
        second.key <- candidates$second[i]
        f <- first.left
        s <- second.left
        f[first.key] <- f[first.key] - 1L
        s[second.key] <- s[second.key] - 1L
        feasible <- all(f <= sum(s) - s) && all(s <= sum(f) - f)
        same.penalty <- if (identical(first.key, second.key)) 1e7 else 0
        feasible.penalty <- if (!feasible) 1e5 else 0
        pair.key <- paste(sort(c(first.key, second.key)), collapse = "\r")
        ordered.key <- paste(first.key, second.key, sep = "\r")
        pair.penalty <- if (identical(first.key, second.key)) {
          0
        } else {
          10 * pair.count[pair.key] + ordered.count[ordered.key]
        }
        same.penalty + feasible.penalty + pair.penalty
      }, numeric(1L))
      chosen <- sample(which(score == min(score)), 1L)
      first.key <- candidates$first[chosen]
      second.key <- candidates$second[chosen]
      out[[b]] <- c(fp[[first.key]][1L], sp[[second.key]][1L])
      key.out[b, ] <- c(first.key, second.key)
      fp[[first.key]] <- fp[[first.key]][-1L]
      sp[[second.key]] <- sp[[second.key]][-1L]
      if (!identical(first.key, second.key)) {
        pair.key <- paste(sort(c(first.key, second.key)), collapse = "\r")
        ordered.key <- paste(first.key, second.key, sep = "\r")
        pair.count[pair.key] <- pair.count[pair.key] + 1L
        ordered.count[ordered.key] <- ordered.count[ordered.key] + 1L
      }
    }

    same <- sum(key.out[, 1L] == key.out[, 2L])
    list(
      block.items = out,
      same = same,
      score = same * 1e6 + sum(pair.count^2) + sum(ordered.count^2)
    )
  }

  best <- NULL
  for (attempt in seq_len(100L)) {
    candidate <- build_once()
    if (is.null(best) || candidate$score < best$score) best <- candidate
    if (candidate$same == 0L) return(candidate$block.items)
  }

  best$block.items
}
#' Compute FCDCM Block Response Probabilities
#' @param zeta Binary statement-level mastery status matrix.
#' @param par Block parameter matrix. A \eqn{B \times 2} matrix with columns
#'   \code{eta0} and \code{etaAB}. Equal-condensation probabilities
#'   (\eqn{\zeta_{A_b}=\zeta_{B_b}}) are fixed at 0.5.
#' @param patterns Block statement-pair matrix. Column 1 is the first
#'   statement \eqn{A_b}; column 2 is the second statement \eqn{B_b}.
#' @return Matrix of probabilities for choosing the first statement in each
#'   block:
#'   \eqn{\eta_{0b}} if \eqn{\zeta_{A_b}<\zeta_{B_b}},
#'   \eqn{0.5} if equal, and
#'   \eqn{0.5+\eta_{ABb}} if \eqn{\zeta_{A_b}>\zeta_{B_b}}.
#' @keywords internal
model.FCDCM <- function(zeta, par, patterns) {
  zeta <- as.matrix(zeta)
  patterns <- as.matrix(patterns)
  storage.mode(zeta) <- "integer"
  storage.mode(patterns) <- "integer"
  par <- fcdcm_normalize_par(par, B = nrow(patterns), arg = "par")
  cpp_model_FCDCM(zeta, par, patterns)
}

get.zeta <- function(alpha, Q.matrix, dcm.type = "DINA") {
  alpha <- as.matrix(alpha)
  Q.matrix <- as.matrix(Q.matrix)
  storage.mode(alpha) <- "numeric"
  storage.mode(Q.matrix) <- "numeric"

  N.person <- nrow(alpha)
  I.states <- nrow(Q.matrix)

  dcm.type <- toupper(dcm.type)
  if (length(dcm.type) != I.states) {
    dcm.type <- rep(dcm.type[1L], I.states)
  }

  eta <- alpha %*% t(Q.matrix)
  required <- rowSums(Q.matrix)
  zeta <- matrix(0L, N.person, I.states)
  is.dina <- dcm.type == "DINA"
  if (any(is.dina)) {
    zeta[, is.dina] <- sweep(eta[, is.dina, drop = FALSE], 2L,
                             required[is.dina], "==")
  }
  if (any(!is.dina)) {
    zeta[, !is.dina] <- eta[, !is.dina, drop = FALSE] > 0
  }

  storage.mode(zeta) <- "integer"
  return(zeta)
}

fcdcm_alpha_profile_logprob <- function(theta, delta1, delta0, alpha.patterns,
                                        eps = 1e-12) {
  theta <- as.numeric(theta)
  alpha.patterns <- as.matrix(alpha.patterns)
  prob.alpha <- distribution_higher_order(theta, delta1, delta0)
  prob.alpha <- pmin(pmax(prob.alpha, eps), 1 - eps)

  log.alpha <- log(prob.alpha) %*% t(alpha.patterns) +
    log1p(-prob.alpha) %*% t(1 - alpha.patterns)
  max.log <- apply(log.alpha, 1L, max)
  log.denom <- max.log +
    log(rowSums(exp(sweep(log.alpha, 1L, max.log, "-"))))
  sweep(log.alpha, 1L, log.denom, "-")
}

fcdcm_alpha_profile_prob <- function(theta, delta1, delta0, alpha.patterns,
                                     eps = 1e-12) {
  exp(fcdcm_alpha_profile_logprob(theta, delta1, delta0, alpha.patterns, eps))
}

fcdcm_latent_support <- function(theta, pi, delta1, delta0, par,
                                 alpha.patterns, zeta.patterns, patterns,
                                 eps = 1e-12) {
  theta <- as.matrix(theta)
  if (ncol(theta) != 1L) {
    theta <- matrix(as.numeric(theta[, 1L]), ncol = 1L)
  }
  pi <- as.vector(pi)
  if (length(pi) != nrow(theta)) {
    stop("'pi' must have one weight per theta support point.", call. = FALSE)
  }
  pi <- pmax(pi, 0)
  pi <- pi / sum(pi)

  support <- cpp_fcdcm_latent_support(
    theta = theta,
    pi = pi,
    delta1 = as.numeric(delta1),
    delta0 = as.numeric(delta0),
    par = par,
    alpha_patterns = alpha.patterns,
    zeta_patterns = zeta.patterns,
    patterns = patterns
  )
  colnames(support$theta) <- colnames(theta)
  colnames(support$alpha) <- colnames(alpha.patterns)
  colnames(support$alpha.weight) <- rownames(alpha.patterns)
  colnames(support$prob) <- rownames(patterns)
  colnames(support$class.prob) <- rownames(patterns)
  support
}

fcdcm_free_par_vec <- function(delta1, delta0, par) {
  par <- fcdcm_normalize_par(par)
  par.free <- par
  par.free[, "eta0"] <- pmin(pmax(par.free[, "eta0"] / 0.5, 1e-8), 1 - 1e-8)
  par.free[, "etaAB"] <- pmin(pmax(par.free[, "etaAB"] / 0.5, 1e-8), 1 - 1e-8)
  c(
    log(pmax(as.numeric(delta1), 1e-8)),
    as.numeric(delta0),
    as.numeric(qlogis(par.free))
  )
}

fcdcm_expand_free_par <- function(par.vec, D, B, dim.names = NULL,
                                  block.names = NULL,
                                  par.names = c("eta0", "etaAB")) {
  idx <- 1L
  delta1 <- exp(par.vec[idx:(idx + D - 1L)])
  idx <- idx + D
  delta0 <- par.vec[idx:(idx + D - 1L)]
  idx <- idx + D
  K <- length(par.names)
  par <- matrix(plogis(par.vec[idx:(idx + B * K - 1L)]),
                nrow = B, ncol = K)
  colnames(par) <- par.names
  if ("eta0" %in% par.names) par[, "eta0"] <- 0.5 * par[, "eta0"]
  if ("etaAB" %in% par.names) par[, "etaAB"] <- 0.5 * par[, "etaAB"]
  par <- fcdcm_normalize_par(par, B = B)
  pars <- list(delta1 = delta1, delta0 = delta0, par = par)
  names(pars$delta1) <- names(pars$delta0) <- dim.names
  rownames(pars$par) <- block.names
  pars
}

fcdcm_normalize_par <- function(par, B = NULL, arg = "par") {
  par <- as.matrix(par)
  storage.mode(par) <- "numeric"
  if (!is.null(B) && nrow(par) != B) {
    stop("'", arg, "' must have ", B, " rows.", call. = FALSE)
  }
  if (ncol(par) != 2L) {
    stop("'", arg, "' must be a B x 2 matrix (eta0, etaAB).",
         call. = FALSE)
  }
  out <- cbind(eta0 = par[, 1L], etaAB = par[, 2L])
  rownames(out) <- rownames(par)
  if (anyNA(out) || any(!is.finite(out))) {
    stop("'", arg, "' contains non-finite values.", call. = FALSE)
  }
  if (any(out[, "eta0"] <= 0 | out[, "eta0"] >= 0.5) ||
      any(out[, "etaAB"] <= 0 | out[, "etaAB"] >= 0.5)) {
    stop("FCDCM parameters require 0 < eta0 < .5 and 0 < etaAB < .5.",
         call. = FALSE)
  }
  out
}
