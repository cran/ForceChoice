################################# fit shared helpers #################################

fc_as_control_list <- function(control, arg) {
  if (is.null(control)) {
    return(list())
  }
  if (!is.list(control)) {
    stop("'", arg, "' must be NULL or a named list.", call. = FALSE)
  }
  if (length(control) > 0L) {
    control_names <- names(control)
    if (is.null(control_names) || any(!nzchar(control_names))) {
      stop("'", arg, "' must be a named list.", call. = FALSE)
    }
  }
  control
}


################################# model name aliases #################################

model_aliases <- list(
  "rasch" = "m1pl", "1pl" = "m1pl", "m1pl" = "m1pl",
  "2pl"  = "m2pl", "m2pl" = "m2pl",
  "3pl"  = "m3pl", "m3pl" = "m3pl",
  "4pl"  = "m4pl", "m4pl" = "m4pl"
)

#' Resolve IRT model type with alias support
#'
#' Accepts both internal codes (\code{"m1pl"}, \code{"m2pl"},
#' \code{"m3pl"}, \code{"m4pl"}) and user-friendly aliases
#' (\code{"Rasch"}, \code{"1PL"}, \code{"2PL"}, \code{"3PL"},
#' \code{"4PL"}).
#'
#' @param model Character scalar; the model identifier.
#' @return Lower-case internal model code.
#' @keywords internal
resolve_model_type <- function(model) {
  if (length(model) != 1L || is.na(model) || !nzchar(as.character(model))) {
    stop("'model' must be a non-missing character scalar.", call. = FALSE)
  }
  model <- tolower(as.character(model))
  resolved <- model_aliases[[model]]
  if (is.null(resolved)) {
    stop(
      "'model' must be one of: ",
      paste(sort(unique(unlist(model_aliases, use.names = FALSE))), collapse = ", "),
      " (or their aliases: Rasch, 1PL, 2PL, 3PL, 4PL).",
      call. = FALSE
    )
  }
  resolved
}

#' Map resolved model name to integer code
#'
#' Converts the canonical model string returned by \code{resolve_model_type()}
#' to the integer code expected by the Stan models (1 = 1PL, 2 = 2PL,
#' 3 = 3PL, 4 = 4PL).
#'
#' @param model Character string (canonical model name).
#' @return Integer model type code.
#' @keywords internal
model_type_to_int <- function(model) {
  switch(model,
    m1pl = 1L, m2pl = 2L, m3pl = 3L, m4pl = 4L,
    stop("Unknown model: ", model, call. = FALSE)
  )
}


################################# auto-detection helpers #################################

#' Resolve the number of latent dimensions
#'
#' Uses an explicit \eqn{D}, or infers it from the number of columns in a
#' supplied Q-matrix. Latent dimensionality is never guessed from the number
#' of items or blocks.
#'
#' @param D User-supplied D or NULL.
#' @param Q.matrix Optional Q-matrix.
#' @param block.items Deprecated internal argument; ignored.
#' @param I Deprecated internal argument; ignored.
#' @return Integer D.
#' @keywords internal
resolve_D <- function(D, Q.matrix = NULL, block.items = NULL, I = NULL) {
  if (!is.null(D)) {
    return(check_integer_scalar(D, "D"))
  }

  if (!is.null(Q.matrix)) {
    Q.matrix <- as.matrix(Q.matrix)
    if (ncol(Q.matrix) >= 1L) {
      return(ncol(Q.matrix))
    }
  }

  stop("'D' must be supplied when 'Q.matrix' is NULL.", call. = FALSE)
}


#' Auto-detect forced-choice response type from data
#'
#' When \code{fc.type} is \code{NULL}, inspects the ranking strings in the
#' data columns to infer whether they represent RANK, MOLE, or PICK.
#'
#' @param fc.type User-supplied fc.type or NULL.
#' @param data The forced-choice data matrix.
#' @param block.items Optional block-item list.
#' @return Character vector of fc.type values (one per block).
#' @keywords internal
resolve_fc_type <- function(fc.type, data, block.items = NULL) {
  data.matrix <- as.matrix(data)
  N.block <- ncol(data.matrix)
  if (!is.null(fc.type)) {
    return(normalize_fc_type(fc.type, N.block))
  }

  if (!is.character(data.matrix)) {
    return(rep("RANK", N.block))
  }

  if (!is.null(block.items) && length(block.items) != N.block) {
    stop("'block.items' must have one element per data column.",
         call. = FALSE)
  }

  inferred <- character(N.block)
  for (b in seq_len(N.block)) {
    x <- trimws(as.character(data.matrix[, b]))
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) {
      stop("Block ", b, " contains no non-missing responses.",
           call. = FALSE)
    }
    widths <- lengths(strsplit(x, ">", fixed = TRUE))
    if (length(unique(widths)) != 1L) {
      stop("Block ", b, " mixes response formats; supply a consistent ",
           "'fc.type' and 'block.items'.", call. = FALSE)
    }
    width <- widths[1L]
    if (is.null(block.items)) {
      parsed <- lapply(strsplit(x, ">", fixed = TRUE), function(parts) {
        suppressWarnings(as.integer(trimws(parts)))
      })
      if (any(vapply(parsed, anyNA, logical(1L)))) {
        stop("Cannot parse item indices in block ", b, ".", call. = FALSE)
      }
      item_sets <- vapply(
        parsed, function(items) paste(sort(unique(items)), collapse = ","),
        character(1L)
      )
      if (length(unique(item_sets)) != 1L ||
          length(unique(parsed[[1L]])) != width) {
        stop("Partial rankings require explicit 'block.items' and ",
             "'fc.type'.", call. = FALSE)
      }
      inferred[b] <- "RANK"
    } else {
      K <- length(block.items[[b]])
      inferred[b] <- if (width == K) {
        "RANK"
      } else if (width == 2L && K > 2L) {
        "MOLE"
      } else if (width == 1L) {
        "PICK"
      } else {
        stop("Block ", b, " response width does not match 'block.items'.",
             call. = FALSE)
      }
    }
  }
  inferred
}


#' Auto-detect block.items from forced-choice ranking data
#'
#' Parses ranking strings like \code{"2>1>3"} to extract the set of item
#' indices in each block.  Only works when data columns contain character
#' ranking strings.
#'
#' @param block.items User-supplied block.items or NULL.
#' @param data The forced-choice data matrix.
#' @return A list of integer vectors (one per block), or NULL if
#'   auto-detection is not possible.
#' @keywords internal
resolve_block_items <- function(block.items, data) {
  data.matrix <- as.matrix(data)
  N.block <- ncol(data.matrix)
  if (!is.null(block.items)) {
    if (!is.list(block.items) || length(block.items) != N.block) {
      stop("'block.items' must be a list with one element per data column.",
           call. = FALSE)
    }
    block.items <- lapply(seq_along(block.items), function(b) {
      items <- block.items[[b]]
      if (!is.numeric(items) || anyNA(items) || any(!is.finite(items)) ||
          any(items != floor(items)) || any(items < 1) ||
          anyDuplicated(items)) {
        stop("'block.items[[", b, "]]' must contain unique positive ",
             "integer indices.", call. = FALSE)
      }
      as.integer(items)
    })
    if (anyDuplicated(unlist(block.items, use.names = FALSE))) {
      stop("Item indices must be unique across 'block.items'.",
           call. = FALSE)
    }
    return(block.items)
  }

  if (!is.character(data.matrix)) return(NULL)

  block.items <- vector("list", N.block)

  for (b in seq_len(N.block)) {
    x <- as.character(data.matrix[, b])
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) {
      stop("Block ", b, " contains no valid data for auto-detection. ",
           "Please supply 'block.items' explicitly.", call. = FALSE)
    }
    parsed <- lapply(strsplit(x, ">", fixed = TRUE), function(parts) {
      suppressWarnings(as.integer(trimws(parts)))
    })
    if (any(vapply(parsed, anyNA, logical(1L)))) {
      stop("Cannot parse item indices from block ", b, " data. ",
           "Please supply 'block.items' explicitly.", call. = FALSE)
    }
    item.sets <- lapply(parsed, function(items) sort(unique(items)))
    keys <- vapply(
      item.sets, function(items) paste(items, collapse = ","), character(1L)
    )
    if (length(unique(keys)) != 1L ||
        any(lengths(item.sets) != lengths(parsed))) {
      stop("Block ", b, " contains partial or inconsistent rankings. ",
           "Please supply 'block.items' explicitly.", call. = FALSE)
    }
    block.items[[b]] <- item.sets[[1L]]
  }

  # Validate: items should form a contiguous set 1..I
  all_items <- sort(unique(unlist(block.items, use.names = FALSE)))
  I <- length(all_items)
  if (!all(all_items == seq_len(I))) {
    stop("Auto-detected block items are not a contiguous set 1..", I, ". ",
         "Please supply 'block.items' explicitly.", call. = FALSE)
  }

  message("block.items auto-detected: ", N.block, " blocks, ",
          I, " unique items")
  block.items
}

fit_common_method_control <- function(control.method) {
  control.method <- fc_as_control_list(control.method, "control.method")

  chains <- as.integer(get_ctrl("chains", 2L, control.method))
  cores <- get_ctrl("cores", chains, control.method)
  cores <- as.integer(cores)
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    stop("'cores' in 'control.method' must be a positive integer.",
         call. = FALSE)
  }

  vis <- get_ctrl("vis", TRUE, control.method)
  if (!is.logical(vis) || length(vis) != 1L || is.na(vis)) {
    stop("'vis' in 'control.method' must be TRUE or FALSE.", call. = FALSE)
  }

  seed <- get_ctrl(
    "seed", sample.int(.Machine$integer.max, 1L), control.method
  )
  if ((!is.numeric(seed) && !is.integer(seed)) ||
      length(seed) != 1L || !is.finite(seed)) {
    stop("'seed' in 'control.method' must be a finite numeric scalar.",
         call. = FALSE)
  }
  seed <- as.integer(seed)
  if (is.na(seed)) {
    stop("'seed' in 'control.method' must be a valid integer.", call. = FALSE)
  }

  list(cores = cores, vis = vis, seed = seed)
}

fit_effective_method_control <- function(control.method, method = list(),
                                         common = list()) {
  control.method <- fc_as_control_list(control.method, "control.method")
  defaults <- utils::modifyList(common, method)
  utils::modifyList(defaults, control.method)
}

fc_stan_method_control <- function(control.method) {
  control.method <- fc_as_control_list(control.method, "control.method")

  iter <- as.integer(get_ctrl("iter", 5000L, control.method))
  stan.method <- list(
    chains = as.integer(get_ctrl("chains", 2L, control.method)),
    iter = iter,
    warmup = as.integer(get_ctrl("warmup", floor(iter / 2), control.method)),
    thin = as.integer(get_ctrl("thin", 1L, control.method)),
    init = get_ctrl("init", "random", control.method),
    algorithm = get_ctrl("algorithm", "HMC", control.method)
  )

  stan.method$algorithm <- match.arg(
    stan.method$algorithm,
    c("NUTS", "HMC", "Fixed_param")
  )

  integer_names <- c("chains", "iter", "warmup", "thin")
  for (nm in integer_names) {
    if (is.na(stan.method[[nm]]) || stan.method[[nm]] < 1L) {
      stop("'", nm, "' in 'control.method' must be a positive integer.",
           call. = FALSE)
    }
  }
  if (stan.method$warmup >= stan.method$iter) {
    stop("'warmup' in 'control.method' must be smaller than 'iter'.",
         call. = FALSE)
  }

  stan.method
}

istem_scalar <- function(value, name, positive = FALSE, arg = "control.model") {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("'", name, "' in '", arg, "' must be a finite numeric scalar.",
         call. = FALSE)
  }
  value <- as.numeric(value)
  if (positive && value <= 0) {
    stop("'", name, "' in '", arg, "' must be positive.", call. = FALSE)
  }
  value
}

istem_clip <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

istem_random_theta <- function(N, D, theta.mu, lower, upper) {
  theta.mu <- rep(as.numeric(theta.mu), length.out = D)
  theta <- matrix(stats::rnorm(N * D, mean = rep(theta.mu, each = N), sd = 1),
                  nrow = N, ncol = D)
  istem_clip(theta, lower, upper)
}

istem_rnorm_bounded <- function(n, mean, sd, lower, upper) {
  istem_clip(stats::rnorm(n, mean = mean, sd = sd), lower, upper)
}

istem_rlnorm_bounded <- function(n, meanlog, sdlog, lower, upper) {
  istem_clip(stats::rlnorm(n, meanlog = meanlog, sdlog = sdlog), lower, upper)
}

istem_random_corr <- function(D, include.corr = TRUE) {
  if (D <= 1L || !isTRUE(include.corr)) {
    return(diag(D))
  }
  z <- matrix(stats::rnorm(D * (D + 1L)), nrow = D)
  corr <- stats::cov2cor(tcrossprod(z))
  diag(corr) <- 1
  corr
}

istem_response_groups <- function(response) {
  response <- as.matrix(response)
  key <- apply(response, 1L, paste, collapse = "\r")
  level <- unique(key)
  group <- match(key, level)
  first <- match(level, key)
  unique.response <- response[first, , drop = FALSE]
  storage.mode(unique.response) <- storage.mode(response)
  list(
    response = unique.response,
    count = as.integer(tabulate(group, nbins = length(level))),
    group = as.integer(group),
    first = as.integer(first),
    G = length(level)
  )
}

istem_response_pair_groups <- function(response, pairs.value) {
  response <- as.matrix(response)
  response.key <- apply(response, 1L, paste, collapse = "\r")
  pair.key <- vapply(pairs.value, function(person) {
    paste(vapply(person, function(block) {
      block <- as.matrix(block)
      paste(nrow(block), ncol(block), paste(as.integer(block), collapse = ","),
            sep = ":")
    }, character(1L)), collapse = "\r")
  }, character(1L))
  key <- paste(response.key, pair.key, sep = "\f")
  level <- unique(key)
  group <- match(key, level)
  first <- match(level, key)
  unique.response <- response[first, , drop = FALSE]
  storage.mode(unique.response) <- storage.mode(response)
  list(
    response = unique.response,
    pairs.value = pairs.value[first],
    count = as.integer(tabulate(group, nbins = length(level))),
    group = as.integer(group),
    first = as.integer(first),
    G = length(level)
  )
}

istem_expand_group_matrix <- function(x, group, row_names = NULL) {
  out <- x[group, , drop = FALSE]
  if (!is.null(row_names)) rownames(out) <- row_names
  out
}

istem_state_weight <- function(state) {
  if (!is.null(state$response.count)) state$response.count else rep(1, nrow(state$theta))
}

istem_default_L <- function(D) {
  switch(as.character(D),
         "1" = 61L, "2" = 31L, "3" = 15L, "4" = 9L, "5" = 7L, 3L)
}

istem_latent_grid_control <- function(control.model, D, L = NULL,
                                      theta.low = NULL, theta.up = NULL) {
  control.model <- fc_as_control_list(control.model, "control.model")
  D <- as.integer(D)
  if (length(D) != 1L || is.na(D) || D < 1L) {
    stop("'D' must be a positive integer.", call. = FALSE)
  }

  if (is.null(L)) {
    L <- get_ctrl("L", NULL, control.model)
  }
  if (is.null(L)) {
    L <- istem_default_L(D)
  }
  L <- as.integer(L)
  if (length(L) != 1L || is.na(L) || !is.finite(L) || L < 2L) {
    stop("'L' in 'control.model' must be an integer of at least 2.",
         call. = FALSE)
  }

  if (is.null(theta.low)) {
    theta.low <- get_ctrl("theta.lower", -6, control.model)
  }
  if (is.null(theta.up)) {
    theta.up <- get_ctrl("theta.upper", 6, control.model)
  }
  theta.low <- as.numeric(theta.low)
  theta.up <- as.numeric(theta.up)
  if (length(theta.low) != 1L || length(theta.up) != 1L ||
      !is.finite(theta.low) || !is.finite(theta.up) ||
      theta.low >= theta.up) {
    stop("'theta.lower' must be smaller than 'theta.upper' in 'control.model'.",
         call. = FALSE)
  }

  list(L = L, theta.lower = theta.low, theta.upper = theta.up)
}

istem_apply_grid_control <- function(method, control.model, D) {
  grid <- istem_latent_grid_control(control.model, D)
  method$L <- grid$L
  method$theta.lower <- grid$theta.lower
  method$theta.upper <- grid$theta.upper
  method
}

istem_effective_grid_model_control <- function(control.model, method) {
  utils::modifyList(
    control.model,
    list(L = method$L,
         theta.lower = method$theta.lower,
         theta.upper = method$theta.upper)
  )
}

istem_method_control <- function(control.method) {
  control.method <- fc_as_control_list(control.method, "control.method")
  list(
    M = as.integer(get_ctrl("M", 10L, control.method)),
    B = as.integer(get_ctrl("B", 20L, control.method)),
    burnin.maxitr = as.integer(get_ctrl("burnin.maxitr", 100L, control.method)),
    maxitr = as.integer(get_ctrl("maxitr", 2000L, control.method)),
    eps1 = get_ctrl("eps1", 1.5, control.method),
    eps2 = get_ctrl("eps2", 0.4, control.method),
    frac1 = get_ctrl("frac1", 0.1, control.method),
    frac2 = get_ctrl("frac2", 0.5, control.method),
    optim.maxit = as.integer(get_ctrl("optim.maxit", 200L, control.method)),
    corr.optim.maxit = as.integer(get_ctrl("corr.optim.maxit", 50L,
                                           control.method)),
    fix.corr = isTRUE(get_ctrl("fix.corr", FALSE, control.method)),
    estimate.se = isTRUE(get_ctrl("estimate.se", TRUE, control.method))
  )
}

istem_check_method_control <- function(method) {
  if (method$M < 2L || method$B < 1L) {
    stop("'M' must be at least 2 and 'B' must be at least 1 for iStEM.",
         call. = FALSE)
  }
  if (method$burnin.maxitr < method$M) {
    stop("'burnin.maxitr' must be at least as large as 'M'.", call. = FALSE)
  }
  if (method$maxitr < method$M) {
    stop("'maxitr' must be at least as large as 'M'.", call. = FALSE)
  }
  if (!is.null(method$L) &&
      (!is.finite(method$L) || as.integer(method$L) < 2L)) {
    stop("'L' in 'control.model' must be an integer of at least 2.",
         call. = FALSE)
  }
  if (!is.null(method$theta.lower) || !is.null(method$theta.upper)) {
    if (!is.finite(method$theta.lower) || !is.finite(method$theta.upper) ||
        method$theta.lower >= method$theta.upper) {
      stop("'theta.lower' must be smaller than 'theta.upper' in 'control.model'.",
           call. = FALSE)
    }
  }
  if (!is.finite(method$corr.optim.maxit) || method$corr.optim.maxit < 1L) {
    stop("'corr.optim.maxit' must be a positive integer.", call. = FALSE)
  }
  invisible(NULL)
}

istem_progress <- function(fmt, ..., width = 120L, newline = FALSE) {
  msg <- sprintf(fmt, ...)
  pad <- strrep(" ", max(0L, as.integer(width) - nchar(msg, type = "width")))
  cat("\r", msg, pad, if (newline) "\n" else "", sep = "")
  flush.console()
  invisible(NULL)
}

istem_attach_corr_cache <- function(C) {
  C <- as.matrix(C)
  chol_sigma <- tryCatch(chol(C), error = function(e) NULL)
  if (is.null(chol_sigma)) {
    return(C)
  }
  attr(C, "chol") <- chol_sigma
  attr(C, "log_diag_sum") <- sum(log(diag(chol_sigma)))
  C
}

istem_safe_corr <- function(x) {
  if (is.matrix(x) && nrow(x) == ncol(x)) {
    C <- x
  } else {
    C <- stats::cor(x)
  }
  C[!is.finite(C)] <- 0
  C <- (C + t(C)) / 2
  diag(C) <- 1
  eig <- tryCatch(eigen(C, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig)) {
    C <- diag(ncol(C))
    return(istem_attach_corr_cache(C))
  }
  eig$values <- pmax(eig$values, 1e-6)
  C <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*%
    t(eig$vectors)
  C <- stats::cov2cor(C)
  C[!is.finite(C)] <- 0
  C <- (C + t(C)) / 2
  diag(C) <- 1
  istem_attach_corr_cache(C)
}

istem_corr_chol <- function(Corr) {
  chol_sigma <- attr(Corr, "chol", exact = TRUE)
  log_diag_sum <- attr(Corr, "log_diag_sum", exact = TRUE)
  if (is.null(chol_sigma) || is.null(log_diag_sum)) {
    Corr <- istem_safe_corr(Corr)
    chol_sigma <- attr(Corr, "chol", exact = TRUE)
    log_diag_sum <- attr(Corr, "log_diag_sum", exact = TRUE)
  }
  list(chol = chol_sigma, log_diag_sum = log_diag_sum)
}

istem_corr_objective <- function(v, theta, theta.mu) {
  D <- ncol(theta)
  C <- diag(D)
  C[lower.tri(C)] <- v
  C <- C + t(C)
  diag(C) <- 1

  chol_sigma <- tryCatch(chol(C), error = function(e) NULL)
  if (is.null(chol_sigma)) {
    return(1e100)
  }

  z <- tryCatch(
    forwardsolve(t(chol_sigma), t(sweep(theta, 2L, theta.mu, "-"))),
    error = function(e) NULL
  )
  if (is.null(z) || any(!is.finite(z))) {
    return(1e100)
  }

  nrow(theta) * 2 * sum(log(diag(chol_sigma))) + sum(z^2)
}

istem_corr_objective_weighted <- function(v, theta, theta.mu, weight) {
  D <- ncol(theta)
  C <- diag(D)
  C[lower.tri(C)] <- v
  C <- C + t(C)
  diag(C) <- 1

  chol_sigma <- tryCatch(chol(C), error = function(e) NULL)
  if (is.null(chol_sigma)) {
    return(1e100)
  }

  z <- tryCatch(
    forwardsolve(t(chol_sigma), t(sweep(theta, 2L, theta.mu, "-"))),
    error = function(e) NULL
  )
  if (is.null(z) || any(!is.finite(z))) {
    return(1e100)
  }

  sum(weight) * 2 * sum(log(diag(chol_sigma))) +
    sum(sweep(z^2, 2L, weight, "*"))
}

istem_update_corr <- function(theta, theta.mu, Corr, maxit) {
  D <- ncol(theta)
  if (D <= 1L) {
    return(matrix(1, 1L, 1L))
  }

  Corr <- istem_safe_corr(Corr)
  x0 <- Corr[lower.tri(Corr)]
  opt <- tryCatch(
    stats::optim(
      par = x0,
      fn = istem_corr_objective,
      theta = theta,
      theta.mu = theta.mu,
      method = "L-BFGS-B",
      lower = rep(-0.995, length(x0)),
      upper = rep(0.995, length(x0)),
      control = list(maxit = maxit)
    ),
    error = function(e) NULL
  )

  if (!is.null(opt) && is.finite(opt$value)) {
    C <- diag(D)
    C[lower.tri(C)] <- opt$par
    C <- C + t(C)
    diag(C) <- 1
    return(istem_safe_corr(C))
  }

  Corr
}

istem_update_corr_weighted <- function(theta, theta.mu, weight, Corr, maxit) {
  D <- ncol(theta)
  if (D <= 1L) {
    return(matrix(1, 1L, 1L))
  }

  Corr <- istem_safe_corr(Corr)
  x0 <- Corr[lower.tri(Corr)]
  opt <- tryCatch(
    stats::optim(
      par = x0,
      fn = istem_corr_objective_weighted,
      theta = theta,
      theta.mu = theta.mu,
      weight = weight,
      method = "L-BFGS-B",
      lower = rep(-0.995, length(x0)),
      upper = rep(0.995, length(x0)),
      control = list(maxit = maxit)
    ),
    error = function(e) NULL
  )

  if (!is.null(opt) && is.finite(opt$value)) {
    C <- diag(D)
    C[lower.tri(C)] <- opt$par
    C <- C + t(C)
    diag(C) <- 1
    return(istem_safe_corr(C))
  }

  Corr
}

istem_batch_var <- function(plist, n) {
  out <- cpp_istem_batch_var(plist, as.integer(n))
  out[!is.finite(out)] <- Inf
  out
}

istem_geweke_z <- function(chain, frac1 = 0.2, frac2 = 0.5) {
  if (is.null(dim(chain))) {
    chain <- matrix(chain, nrow = 1L)
  }
  z <- coda::geweke.diag(
    coda::mcmc(t(chain)),
    frac1 = frac1,
    frac2 = frac2
  )$z
  z[!is.finite(z)] <- 0
  z
}

istem_kernel <- function(state, B, theta_grid_length, theta_lower, theta_upper,
                         fix.corr, corr_optim_maxit, update_parameters, param_vec,
                         sample_theta, logLik_fun = NULL) {
  state$Corr <- istem_attach_corr_cache(state$Corr)
  N <- nrow(state$theta)
  D <- ncol(state$theta)
  npar <- length(param_vec(state))
  pmatrix <- matrix(NA_real_, npar, B)
  logLik.trace <- rep(NA_real_, B)
  theta_sum <- matrix(0, N, D)
  theta_sq_sum <- matrix(0, N, D)

  for (b in seq_len(B)) {
    sampled <- sample_theta(
      state = state,
      theta_grid_length = theta_grid_length,
      theta_lower = theta_lower,
      theta_upper = theta_upper
    )
    state$theta <- sampled$theta

    if (!fix.corr && D > 1L) {
      weight <- istem_state_weight(state)
      state$Corr <- if (is.null(state$response.count)) {
        istem_update_corr(
          theta = state$theta,
          theta.mu = state$theta_mu,
          Corr = state$Corr,
          maxit = corr_optim_maxit
        )
      } else {
        istem_update_corr_weighted(
          theta = state$theta,
          theta.mu = state$theta_mu,
          weight = weight,
          Corr = state$Corr,
          maxit = corr_optim_maxit
        )
      }
    }

    state <- update_parameters(state)
    pmatrix[, b] <- param_vec(state)
    if (!is.null(logLik_fun)) {
      logLik.trace[b] <- logLik_fun(state)
    }
    theta_sum <- theta_sum + state$theta
    theta_sq_sum <- theta_sq_sum + state$theta^2
  }

  list(
    pmatrix = pmatrix,
    state = state,
    logLik.trace = logLik.trace,
    theta_sum = theta_sum,
    theta_sq_sum = theta_sq_sum
  )
}

istem_trace_logLik_values <- function(x) {
  value <- x$logLik.trace %||% x$logLik %||% x$log_lik %||% x$loglik
  if (is.null(value)) return(NA_real_)
  if (is.data.frame(value)) {
    value <- value$logLik %||% value$log_lik %||% value$loglik
  }
  value <- as.numeric(value)
  if (length(value) == 0L) NA_real_ else value
}

istem_run_batches <- function(state, method, N, vis, label, run_batch) {
  M <- method$M
  B <- method$B
  burnin.maxitr <- method$burnin.maxitr
  maxitr <- method$maxitr
  eps1 <- method$eps1
  eps2 <- method$eps2
  frac1 <- method$frac1
  frac2 <- method$frac2

  plist <- vector("list", M)
  trace.logLik <- numeric(0L)
  batch.extra.names <- NULL
  batch.extra <- list()

  if (vis) {
    cat("\n", label, " iStEM burn-in phase:\n", sep = "")
  }

  for (m in seq_len(M)) {
    if (vis) {
      istem_progress("  burn-in batch %04d/%04d", m, M)
    }
    x <- run_batch(state)
    state <- x$state
    plist[[m]] <- x$pmatrix
    trace.logLik <- c(trace.logLik, istem_trace_logLik_values(x))
    if (is.null(batch.extra.names)) {
      batch.extra.names <- setdiff(names(x), c("state", "pmatrix", "logLik",
                                               "log_lik", "loglik",
                                               "logLik.trace"))
      batch.extra <- stats::setNames(vector("list", length(batch.extra.names)),
                                     batch.extra.names)
      for (nm in batch.extra.names) batch.extra[[nm]] <- vector("list", M)
    }
    for (nm in batch.extra.names) batch.extra[[nm]][[m]] <- x[[nm]]
  }

  z <- istem_geweke_z(Reduce(cbind, plist), frac1 = frac1, frac2 = frac2)
  d.hat <- istem_batch_var(plist, n = M)

  burn.in.size <- 0L
  total.number.of.batch <- M
  m <- M
  while (sum(z^2) >= length(z) * eps1 &&
         m < burnin.maxitr) {

    x <- run_batch(state)
    state <- x$state
    trace.logLik <- c(trace.logLik, istem_trace_logLik_values(x))
    plist <- c(plist[-1L], list(x$pmatrix))
    for (nm in batch.extra.names) {
      batch.extra[[nm]] <- c(batch.extra[[nm]][-1L], list(x[[nm]]))
    }

    z <- istem_geweke_z(Reduce(cbind, plist), frac1 = frac1, frac2 = frac2)
    d.hat <- istem_batch_var(plist, n = M)
    burn.in.size <- burn.in.size + B
    total.number.of.batch <- total.number.of.batch + 1L
    m <- m + 1L

    if (vis) {
      istem_progress(
        "  burn-in batch %04d | Geweke mean z^2 %10.4f <= %-10.4f",
        m,
        sum(z^2) / length(z),
        eps1
      )
    }
  }

  burnin.z2.per.par <- sum(z^2) / length(z)
  burnin.max.Nu <- max(d.hat * N, na.rm = TRUE)
  burnin.converged <- burnin.z2.per.par < eps1
  if (!burnin.converged && m >= burnin.maxitr) {
    warning(
      label,
      " iStEM burn-in stopped at 'burnin.maxitr' before the Geweke criterion was met: ",
      "Geweke mean z^2 = ", signif(burnin.z2.per.par, 5),
      ", criterion = ", signif(eps1, 5), ".",
      call. = FALSE
    )
  }

  stable.plist <- plist
  stable.extra <- batch.extra
  n <- M
  d.hat <- istem_batch_var(stable.plist, n = n)

  if (vis) {
    cat("\nAfter burn-in phase:\n")
  }

  while (max(d.hat * N, na.rm = TRUE) >= eps2 && n < maxitr) {
    if (vis) {
      istem_progress(
        "  valid batch %04d | max N*u %10.4f <= %-10.4f",
        n,
        max(d.hat * N, na.rm = TRUE),
        eps2
      )
    }

    x <- run_batch(state)
    state <- x$state
    trace.logLik <- c(trace.logLik, istem_trace_logLik_values(x))
    stable.plist[[n + 1L]] <- x$pmatrix
    for (nm in batch.extra.names) stable.extra[[nm]][[n + 1L]] <- x[[nm]]
    n <- n + 1L
    total.number.of.batch <- total.number.of.batch + 1L
    d.hat <- istem_batch_var(stable.plist, n = n)
  }

  if (vis) {
    istem_progress(
      "  valid batch %04d | max N*u %10.4f <= %-10.4f",
      n,
      max(d.hat * N, na.rm = TRUE),
      eps2,
      newline = TRUE
    )
    cat("Length of final MC chain = ", n * B, "\n", sep = "")
  }

  final.max.Nu <- max(d.hat * N, na.rm = TRUE)
  final.converged <- final.max.Nu < eps2
  if (!final.converged && n >= maxitr) {
    warning(
      label,
      " iStEM stopped at 'maxitr' before the Monte Carlo error criterion was met: ",
      "max N*u = ", signif(final.max.Nu, 5),
      ", criterion = ", signif(eps2, 5), ".",
      call. = FALSE
    )
  }

  chain.length <- n * B
  stable.start <- max(1L, length(trace.logLik) - M * B + 1L)
  stable.start.batch <- ceiling(stable.start / B)
  logLik.trace <- data.frame(
    iter = seq_along(trace.logLik),
    logLik = trace.logLik
  )
  list(
    state = state,
    stable.plist = stable.plist,
    stable.extra = stable.extra,
    stable.start = stable.start,
    stable.start.batch = stable.start.batch,
    chain.length = chain.length,
    n.batch = n,
    iStEM = list(
      M = M,
      B = B,
      burn.in.size = burn.in.size,
      total.number.of.batch = total.number.of.batch,
      final.chain = chain.length,
      burnin.converged = burnin.converged,
      converged = final.converged,
      burnin.z2.per.par = burnin.z2.per.par,
      burnin.max.Nu = burnin.max.Nu,
      final.max.Nu = final.max.Nu,
      geweke.z = z,
      batch.var = d.hat,
      stable.start = stable.start,
      stable.start.batch = stable.start.batch,
      logLik.trace = logLik.trace,
      itemlist = stable.plist
    )
  )
}

istem_run <- function(state, method, N, vis, label,
                      update_parameters, param_vec,
                      sample_theta, logLik_fun = NULL) {
  L <- method$L
  theta_lower <- method$theta.lower
  theta_upper <- method$theta.upper
  corr_optim_maxit <- method$corr.optim.maxit
  fix.corr <- method$fix.corr
  estimate.se <- method$estimate.se

  run <- istem_run_batches(
    state = state,
    method = method,
    N = N,
    vis = vis,
    label = label,
    run_batch = function(state) {
      x <- istem_kernel(
        state = state,
        B = method$B,
        theta_grid_length = L,
        theta_lower = theta_lower,
        theta_upper = theta_upper,
        fix.corr = fix.corr,
        corr_optim_maxit = corr_optim_maxit,
        update_parameters = update_parameters,
        param_vec = param_vec,
        sample_theta = sample_theta,
        logLik_fun = logLik_fun
      )
      x
    }
  )

  chain_mat <- Reduce(cbind, run$stable.plist)
  chain_mean <- rowMeans(chain_mat)
  chain_sd <- if (estimate.se && ncol(chain_mat) > 1L) {
    apply(chain_mat, 1L, stats::sd)
  } else {
    rep(NA_real_, length(chain_mean))
  }
  chain_rhat <- rep(NA_real_, length(chain_mean))  # iStEM is not MCMC; Rhat not applicable

  theta_sum <- Reduce("+", run$stable.extra$theta_sum)
  theta_sq_sum <- Reduce("+", run$stable.extra$theta_sq_sum)
  chain_length <- run$chain.length
  theta_est <- theta_sum / chain_length
  theta_second <- theta_sq_sum / chain_length
  theta_se <- sqrt(pmax(theta_second - theta_est^2, 0))

  list(
    state = run$state,
    chain.mean = chain_mean,
    chain.sd = chain_sd,
    chain.rhat = chain_rhat,
    chain.mat = chain_mat,
    theta.est = theta_est,
    theta.se = theta_se,
    iStEM = utils::modifyList(
      run$iStEM,
      list(
      accept.rate = NA_real_,
        L = L
      )
    )
  )
}

istem_loglik_grid <- function(state, L, theta.lower, theta.upper) {
  D <- ncol(state$theta)
  make_quadrature_grid_mvn(
    state$Corr, state$theta_mu, D, L, theta.lower, theta.upper
  )
}

istem_loglik_binary_sum <- function(prob, response, count, pi) {
  response <- as.matrix(response)
  lik <- cpp_loglik_binary(
    prob = prob,
    response = matrix(as.integer(response), nrow = nrow(response),
                      ncol = ncol(response)),
    pi = as.vector(pi)
  )
  loglik_apply_group_count(lik, count)$logLik
}

istem_loglik_indexed_sum <- function(prob, response, count, pi, block_sizes,
                                     response_base) {
  response <- as.matrix(response)
  lik <- cpp_loglik_indexed(
    prob = prob,
    response = matrix(as.integer(response), nrow = nrow(response),
                      ncol = ncol(response)),
    pi = as.vector(pi),
    block_sizes = as.integer(block_sizes),
    response_base = as.integer(response_base)
  )
  loglik_apply_group_count(lik, count)$logLik
}

mirt_istem_loglik_trace <- function(state, L, theta.lower, theta.upper) {
  grid <- istem_loglik_grid(state, L, theta.lower, theta.upper)
  prob <- model.MIRT(grid$theta.norm, state$par)
  istem_loglik_binary_sum(prob, state$response, state$response.count, grid$pi)
}

mgpcm_istem_loglik_trace <- function(state, L, theta.lower, theta.upper) {
  grid <- istem_loglik_grid(state, L, theta.lower, theta.upper)
  prob <- model.MGPCM(grid$theta.norm, state$par)
  istem_loglik_indexed_sum(
    prob, state$response, state$response.count, grid$pi,
    block_sizes = state$length.poly, response_base = 0L
  )
}

mggum_istem_loglik_trace <- function(state, L, theta.lower, theta.upper) {
  grid <- istem_loglik_grid(state, L, theta.lower, theta.upper)
  prob <- model.MGGUM(grid$theta.norm, state$par)
  istem_loglik_indexed_sum(
    prob, state$response, state$response.count, grid$pi,
    block_sizes = state$length.poly, response_base = 0L
  )
}

fcmirt_istem_loglik_trace <- function(state, L, theta.lower, theta.upper) {
  grid <- istem_loglik_grid(state, L, theta.lower, theta.upper)
  prob <- model.FCMIRT(grid$theta.norm, state$par,
                       state$patterns.total, state$patterns)
  istem_loglik_indexed_sum(
    prob, state$response, state$response.count, grid$pi,
    block_sizes = vapply(state$patterns, nrow, integer(1L)),
    response_base = 1L
  )
}

fcggum_istem_loglik_trace <- function(state, L, theta.lower, theta.upper) {
  grid <- istem_loglik_grid(state, L, theta.lower, theta.upper)
  prob <- model.FCGGUM(grid$theta.norm, state$par,
                       state$patterns.total, state$patterns)
  istem_loglik_indexed_sum(
    prob, state$response, state$response.count, grid$pi,
    block_sizes = vapply(state$patterns, nrow, integer(1L)),
    response_base = 1L
  )
}

tirt_istem_loglik_trace <- function(state, L, theta.lower, theta.upper) {
  grid <- istem_loglik_grid(state, L, theta.lower, theta.upper)
  gamma.matrix <- tirt_istem_gamma_matrix(
    state$gamma, rep(NA_real_, length(state$gamma)),
    state$pairs.matrix, state$I.states
  )$est
  prob <- model.TIRT(grid$theta.norm, state$lambda, state$psi^2,
                     gamma.matrix, state$Q.matrix, state$pairs.matrix)

  response <- as.matrix(state$response.unique)
  block.pair.counts <- vapply(seq_along(state$block.items), function(b) {
    k <- length(state$block.items[[b]])
    as.integer(switch(state$fc.type[b],
      RANK = choose(k, 2),
      MOLE = 2L * k - 3L,
      PICK = k - 1L,
      stop("'fc.type' must be RANK, MOLE, or PICK.", call. = FALSE)
    ))
  }, integer(1L))

  if (!is.null(state$pairs.value.unique) && !all(state$fc.type == "RANK")) {
    lik <- cpp_loglik_binary_person_pairs(
      prob = prob,
      response = matrix(as.integer(response), nrow = nrow(response),
                        ncol = ncol(response)),
      pi = as.vector(grid$pi),
      pairs_matrix = state$pairs.matrix,
      pairs_value = state$pairs.value.unique,
      block_pair_counts = block.pair.counts
    )
  } else {
    lik <- cpp_loglik_binary(
      prob = prob,
      response = matrix(as.integer(response), nrow = nrow(response),
                        ncol = ncol(response)),
      pi = as.vector(grid$pi)
    )
  }
  loglik_apply_group_count(lik, state$response.count)$logLik
}

istem_prepare_01_q <- function(I, D, Q.matrix, triangular = TRUE,
                               require.row = FALSE) {
  if (!is.null(Q.matrix)) {
    Q.matrix <- as.matrix(Q.matrix)
    if (nrow(Q.matrix) != I || ncol(Q.matrix) != D ||
        anyNA(Q.matrix) || !all(Q.matrix %in% c(0, 1))) {
      stop("'Q.matrix' must be an I x D matrix containing only 0 and 1.",
           call. = FALSE)
    }
  }

  if (is.null(Q.matrix)) {
    if (triangular && D > I) {
      stop("Triangular identification requires 'D' <= number of items.",
           call. = FALSE)
    }
    Q.matrix <- matrix(1, I, D)
    if (triangular && D > 1L) {
      for (k in seq_len(D)) {
        r <- I - D + k
        for (j in seq_len(D)) {
          if (k > D - j + 1L) {
            Q.matrix[r, j] <- 0
          }
        }
      }
    }
  }

  if (require.row && any(rowSums(Q.matrix) < 1L)) {
    stop("Each item must measure at least one trait (each row of 'Q.matrix' needs a 1).",
         call. = FALSE)
  }

  storage.mode(Q.matrix) <- "numeric"
  Q.matrix
}

istem_prepare_signed_q <- function(I, D, Q.matrix, require.row = FALSE) {
  if (!is.null(Q.matrix)) {
    Q.matrix <- as.matrix(Q.matrix)
    if (nrow(Q.matrix) != I || ncol(Q.matrix) != D ||
        anyNA(Q.matrix) || !all(Q.matrix %in% c(-1, 0, 1))) {
      stop("'Q.matrix' must be an I x D matrix containing only -1, 0, and 1.",
           call. = FALSE)
    }
  }
  if (is.null(Q.matrix)) {
    Q.matrix <- matrix(1, I, D)
  }
  if (require.row && any(rowSums(abs(Q.matrix)) < 1L)) {
    stop("Each item must measure at least one trait (each row of 'Q.matrix' needs a non-zero entry).",
         call. = FALSE)
  }
  storage.mode(Q.matrix) <- "numeric"
  Q.matrix
}

istem_prepare_response <- function(response, binary = FALSE) {
  response <- as.matrix(response)
  if (!is.numeric(response) || length(dim(response)) != 2L ||
      nrow(response) < 1L || ncol(response) < 1L) {
    stop("'data' must be a non-empty numeric matrix.", call. = FALSE)
  }
  if (anyNA(response) || any(!is.finite(response))) {
    stop("'data' must not contain missing or non-finite responses.",
         call. = FALSE)
  }
  if (any(response != floor(response))) {
    stop("'data' responses must be integer-valued.", call. = FALSE)
  }
  if (binary) {
    if (any(!response %in% c(0, 1))) {
      stop("'data' must contain only binary responses coded 0 and 1.",
           call. = FALSE)
    }
  } else if (any(response < 0)) {
    stop("'data' responses must be non-negative integers.", call. = FALSE)
  }
  storage.mode(response) <- "integer"
  response
}

istem_prepare_length_poly <- function(response, length.poly = NULL) {
  response <- istem_prepare_response(response, binary = FALSE)
  observed.length.poly <- apply(response, 2, max, na.rm = TRUE) + 1L
  I <- ncol(response)
  if (is.null(length.poly)) {
    length.poly <- observed.length.poly
  } else {
    if (length(length.poly) != I) {
      if (length(length.poly) == 1L) {
        length.poly <- rep(length.poly, I)
      } else {
        stop("'length.poly' must be a scalar or have length equal to ncol(response).",
             call. = FALSE)
      }
    }
    if (!is.numeric(length.poly) || anyNA(length.poly) ||
        any(!is.finite(length.poly)) || any(length.poly != floor(length.poly))) {
      stop("'length.poly' must contain finite integers.", call. = FALSE)
    }
    length.poly <- as.integer(length.poly)
    if (any(length.poly < 2L)) {
      stop("'length.poly' must contain integers >= 2.", call. = FALSE)
    }
    if (any(observed.length.poly > length.poly)) {
      stop("'length.poly' is smaller than the largest observed response category for at least one item.",
           call. = FALSE)
    }
  }
  if (any(length.poly < 2L)) {
    stop("Every item must have at least two response categories.",
         call. = FALSE)
  }
  as.integer(length.poly)
}

istem_named_theta <- function(theta, theta.se, response, D) {
  colnames(theta) <- colnames(theta.se) <- paste0("Dim.", seq_len(D))
  rownames(theta) <- rownames(theta.se) <- rownames(response)
  theta.Rhat <- matrix(NA_real_, nrow(theta), D, dimnames = dimnames(theta))
  list(est = theta, se = theta.se, Rhat = theta.Rhat)
}

istem_param_corr_vec <- function(Corr, include_corr) {
  if (include_corr) Corr[lower.tri(Corr)] else numeric(0)
}

istem_unpack_corr <- function(pv, pv_se, idx, D, Corr, include_corr) {
  Corr.se <- matrix(NA_real_, D, D)
  if (D == 1L) {
    Corr <- matrix(1, 1, 1)
    Corr.se <- matrix(0, 1, 1)
  } else if (include_corr) {
    n_corr <- D * (D - 1L) / 2L
    corr_vals <- pv[idx:(idx + n_corr - 1L)]
    corr_se <- pv_se[idx:(idx + n_corr - 1L)]
    idx <- idx + n_corr
    Corr <- diag(D)
    Corr[lower.tri(Corr)] <- corr_vals
    Corr <- t(Corr)
    Corr[lower.tri(Corr)] <- corr_vals
    Corr.se[lower.tri(Corr.se)] <- corr_se
    Corr.se <- t(Corr.se)
    Corr.se[lower.tri(Corr.se)] <- corr_se
    diag(Corr.se) <- 0
    Corr <- istem_safe_corr(Corr)
  } else {
    Corr.se[,] <- NA_real_
    diag(Corr.se) <- 0
  }
  list(Corr = Corr, Corr.se = Corr.se, idx = idx)
}

istem_unpack_corr_rhat <- function(pv_rhat, idx, D, include_corr) {
  if (is.null(pv_rhat)) return(matrix(NA_real_, D, D))
  Corr.rhat <- matrix(NA_real_, D, D)
  diag(Corr.rhat) <- NA_real_
  if (D > 1L && include_corr) {
    n_corr <- D * (D - 1L) / 2L
    corr_rhat <- pv_rhat[idx:(idx + n_corr - 1L)]
    Corr.rhat[lower.tri(Corr.rhat)] <- corr_rhat
    Corr.rhat <- t(Corr.rhat)
    Corr.rhat[lower.tri(Corr.rhat)] <- corr_rhat
  }
  Corr.rhat
}


################################# Forced-choice iStEM helpers #################################

istem_prepare_fc_data <- function(data, block.items = NULL, fc.type = "RANK") {
  fc.type.valid <- c("RANK", "MOLE", "PICK")
  fc.type <- toupper(fc.type)
  if (!all(fc.type %in% fc.type.valid)) {
    stop("'fc.type' must be one of: ", paste(fc.type.valid, collapse = ", "),
         call. = FALSE)
  }

  data.matrix <- as.matrix(data)
  N.block.data <- ncol(data.matrix)
  if (is.null(N.block.data) || N.block.data < 1L) {
    stop("'data' must have at least one column/block.", call. = FALSE)
  }
  if (length(fc.type) != N.block.data) {
    fc.type <- rep(fc.type[1L], N.block.data)
  }

  is.pattern.index.data <- is.numeric(data.matrix) || is.integer(data.matrix)
  if (is.null(block.items)) {
    if (is.pattern.index.data) {
      stop("'block.items' is required when 'data' is an integer matrix of pattern indices.",
           call. = FALSE)
    }
    if (any(fc.type != "RANK")) {
      stop("'block.items' is required for MOLE/PICK data because partial rankings do not identify all block items.",
           call. = FALSE)
    }
    block.items <- get.block.items.from.data(data)
  }

  if (is.pattern.index.data) {
    response <- matrix(as.integer(data.matrix),
                       nrow = nrow(data.matrix),
                       ncol = ncol(data.matrix),
                       dimnames = dimnames(data.matrix))
  } else {
    response <- get.response.from.data(data, block.items, fc.type = fc.type)
  }

  if (!is.list(block.items)) {
    stop("'block.items' must be a list.", call. = FALSE)
  }
  N.block <- length(block.items)
  all_items <- unlist(block.items, use.names = FALSE)
  if (anyDuplicated(all_items)) {
    stop("'block.items' contains duplicate items.", call. = FALSE)
  }
  I <- length(all_items)
  if (!setequal(all_items, seq_len(I))) {
    stop("'block.items' must contain exactly the item indices 1 through ", I,
         ".", call. = FALSE)
  }
  if (any(vapply(block.items, length, integer(1L)) < 2L)) {
    stop("Each block in 'block.items' must contain at least two items.",
         call. = FALSE)
  }

  response <- as.matrix(response)
  if (ncol(response) != N.block) {
    stop("'response' must have ", N.block, " columns (one per block).",
         call. = FALSE)
  }
  if (length(fc.type) != N.block) {
    fc.type <- rep(fc.type[1L], N.block)
  }

  patterns.total <- vector("list", N.block)
  patterns <- vector("list", N.block)
  for (b in seq_len(N.block)) {
    patterns.total[[b]] <- get_permutations(block.items[[b]],
                                            length(block.items[[b]]))
    patterns[[b]] <- get_permutations(
      block.items[[b]],
      switch(fc.type[b],
             RANK = length(block.items[[b]]),
             MOLE = 2L,
             PICK = 1L)
    )
  }

  for (b in seq_len(N.block)) {
    if (anyNA(response[, b]) ||
        any(response[, b] < 1L | response[, b] > nrow(patterns[[b]]))) {
      stop("'data' contains invalid pattern indices for block ", b,
           "; valid values are 1 through ", nrow(patterns[[b]]), ".",
           call. = FALSE)
    }
  }

  storage.mode(response) <- "integer"

  list(
    response = response,
    block.items = block.items,
    fc.type = fc.type,
    patterns = patterns,
    patterns.total = patterns.total,
    N = nrow(response),
    B = N.block,
    I = I,
    all.items = all_items
  )
}
