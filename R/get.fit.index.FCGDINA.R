#' @describeIn get.fit.index FCGDINA model: forced-choice CDM with nominal
#'   binary expansion; within-block pairs excluded.
#' @export
get.fit.index.FCGDINA <- function(object, ...) {
  response <- object$response
  block.items <- object$block.items
  patterns <- object$patterns

  prep <- fcgdina_prepare(
    data = response,
    Q.matrix = object$Q.matrix,
    model = object$model,
    block.items = block.items,
    fc.type = object$fc.type
  )

  N <- nrow(response)
  B <- length(block.items)
  npar <- object$npar
  logLik.obj <- object$logLik
  pi <- as.vector(attr(logLik.obj, "pi"))

  n_cat <- vapply(patterns, nrow, integer(1L))
  n_bin <- sum(n_cat - 1L)
  response.bin <- matrix(0L, N, n_bin)
  indicator.block <- integer(n_bin)

  idx.resp <- 0L
  for (b in seq_len(B)) {
    Kb <- n_cat[b]
    if (Kb > 1L) {
      for (j in seq_len(Kb - 1L)) {
        response.bin[, idx.resp + j] <- as.integer(response[, b] == j)
        indicator.block[idx.resp + j] <- b
      }
      idx.resp <- idx.resp + (Kb - 1L)
    }
  }

  par.vec <- unlist(object$delta$est, use.names = FALSE)

  expand.free.par <- function(par.vec, ...) {
    delta.lst <- vector("list", prep$I.states)
    off <- 0L
    for (i in seq_len(prep$I.states)) {
      nc <- prep$delta.len[i]
      delta.lst[[i]] <- par.vec[(off + 1L):(off + nc)]
      off <- off + nc
    }
    delta.lst
  }

  prob.full.fun <- function(par.vec, ...) {
    delta.lst <- expand.free.par(par.vec)
    fcgdina_prob_class_link(prep, delta.lst, link = fcgdina_delta_link(object))
  }

  loglik.fun <- function(par.vec, ...) {
    prob.full <- prob.full.fun(par.vec)
    cpp_fcgdina_marginal_loglik(
      prob.full, prep$response, prep$patterns, pi,
      prep$N, prep$B, prep$C)
  }

  prob.fun <- function(par.vec, ...) {
    prob.full <- prob.full.fun(par.vec)
    prob.out <- matrix(0, nrow(prob.full), n_bin)

    idx.full <- 0L
    idx.bin <- 0L
    for (b in seq_len(B)) {
      Kb <- n_cat[b]
      if (Kb > 1L) {
        prob.out[, (idx.bin + 1L):(idx.bin + Kb - 1L)] <-
          prob.full[, (idx.full + 1L):(idx.full + Kb - 1L)]
        idx.bin <- idx.bin + (Kb - 1L)
      }
      idx.full <- idx.full + Kb
    }
    prob.out
  }

  good.of.fit(
    par.vec = par.vec,
    loglik.fun = loglik.fun,
    prob.fun = prob.fun,
    response = response.bin,
    npar = npar,
    pi = pi,
    response.type = "nominal",
    nominal.groups = indicator.block,
    ...
  )
}
