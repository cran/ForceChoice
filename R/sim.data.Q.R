
sim.data.Q.CDM <- function(D, I.states, single=TRUE){
  if(D < 2)
    stop("The number of attributes (D) must be more than 1.\n", call. = FALSE)
  if(I.states < 2)
    stop("The number of items (I.states) must be more than 1.\n", call. = FALSE)

  alpha <- attributepattern(D)
  D.sum <- rowSums(alpha)
  if(single){
    idx <- which(D.sum <= 1 & D.sum > 0)
    # Distribute items evenly across attributes to avoid unbalanced Q-matrices
    # that produce many uninformative within-attribute blocks.
    n.attr <- length(idx)
    base <- I.states %/% n.attr
    remainder <- I.states %% n.attr
    n.per.attr <- rep(base, n.attr)
    if (remainder > 0L) n.per.attr[seq_len(remainder)] <- n.per.attr[seq_len(remainder)] + 1L
    Q <- alpha[rep(idx, times = n.per.attr), , drop = FALSE]
    Q <- Q[sample(nrow(Q)), , drop = FALSE]  # shuffle rows
  }else{
    Q <- alpha[sample(which(D.sum <= 3 & D.sum > 0), I.states, replace = TRUE), ]
  }

  rownames(Q) <- paste0("S", 1:I.states)
  return(Q)
}

attributepattern <- function(D) {
  patterns <- lapply(0:D, function(m) {
    if (m == 0) {
      matrix(0L, 1, D)
    } else {
      t(apply(combn(D, m), 2, function(idx) {
        x <- integer(D)
        x[idx] <- 1L
        x
      }))
    }
  })

  out <- do.call(rbind, patterns)
  colnames(out) <- paste0("A", seq_len(D))
  out
}

sim.Q.matrix.FC <- function(I.states, D, block.items, allow.negative = FALSE) {
  if (length(I.states) != 1L || is.na(I.states) || I.states < 1L) {
    stop("'I.states' must be a positive integer.", call. = FALSE)
  }
  if (length(D) != 1L || is.na(D) || D < 1L) {
    stop("'D' must be a positive integer.", call. = FALSE)
  }
  if (length(allow.negative) != 1L || is.na(allow.negative)) {
    stop("'allow.negative' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(block.items) || length(block.items) < 1L) {
    stop("'block.items' must be a non-empty list.", call. = FALSE)
  }

  I.states <- as.integer(I.states)
  D <- as.integer(D)

  block.items <- lapply(seq_along(block.items), function(b) {
    items.b <- block.items[[b]]
    if (!is.numeric(items.b) && !is.integer(items.b)) {
      stop("Each element of 'block.items' must contain item indices.", call. = FALSE)
    }
    items.b.integer <- as.integer(items.b)
    if (length(items.b.integer) < 1L || anyNA(items.b.integer) ||
        any(items.b != items.b.integer) ||
        any(items.b.integer < 1L | items.b.integer > I.states)) {
      stop("'block.items' must contain valid item indices from 1 to 'I.states'.", call. = FALSE)
    }
    items.b.integer
  })

  all.items <- unlist(block.items, use.names = FALSE)
  if (anyDuplicated(all.items)) {
    stop("'block.items' contains duplicate items.", call. = FALSE)
  }
  if (!setequal(all.items, seq_len(I.states))) {
    stop("'block.items' must contain exactly all item indices from 1 to 'I.states'.", call. = FALSE)
  }

  Q.matrix <- matrix(0, nrow = I.states, ncol = D)
  for (b in seq_along(block.items)) {
    items.b <- block.items[[b]]
    block.traits <- sample(seq_len(D), length(items.b),
                           replace = length(items.b) > D)
    q.values <- if (isTRUE(allow.negative)) {
      sample(c(-1, 1), length(items.b), replace = TRUE)
    } else {
      rep(1, length(items.b))
    }
    Q.matrix[cbind(items.b, block.traits)] <- q.values
  }

  storage.mode(Q.matrix) <- "numeric"
  Q.matrix
}
