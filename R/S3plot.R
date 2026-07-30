#' @title S3 Methods: plot
#'
#' @description
#' Draws log-likelihood trace plots for fitted \pkg{ForceChoice} model
#' objects. Stan, iStEM, and EM fits are visualised through the same
#' convergence target: the total log-likelihood recorded across iterations
#' of the corresponding estimation algorithm.
#'
#' @param x A fitted model object of class \code{"MIRT"},
#'   \code{"MGPCM"}, \code{"MGGUM"}, \code{"FCMIRT"}, \code{"FCDCM"},
#'   \code{"FCGGUM"}, \code{"FCGDINA"}, or \code{"TIRT"}.
#' @param y Ignored.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{NULL}. Plots are produced as a side
#'   effect on the current graphics device.
#'
#' @details
#' The horizontal axis is the saved Stan iteration, iStEM iteration, or EM
#' iteration. The vertical axis is the summed log-likelihood. Stan plots use
#' \code{log_lik} or \code{log_lik_group} from generated quantities and draw
#' all chains in one panel. iStEM plots use the iteration-level
#' \code{iStEM$logLik.trace}, and the stable chain is the final
#' \code{M * B} iStEM iterations. EM plots use \code{EM$logLik.trace}. For
#' Stan and iStEM fits, a vertical dashed line marks the stable-chain
#' boundary. EM fits do not have a stable-chain boundary.
#'
#' @name plot
NULL

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Internal helpers ----

.fc_trace_palette <- function(n) {
  rep_len(c("#4F7FA7", "#B07A7A", "#6F9E7A", "#9A8F68",
            "#7B6FA6", "#6B9CA0"), n)
}

.fc_finite_ylim <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    stop("No finite log-likelihood values are available for plotting.",
         call. = FALSE)
  }
  r <- range(x)
  if (diff(r) == 0) r <- r + c(-0.5, 0.5)
  pad <- diff(r) * 0.05
  r + c(-pad, pad)
}

.fc_finite_xlim <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    stop("No finite iteration values are available for plotting.",
         call. = FALSE)
  }
  r <- range(x)
  if (diff(r) == 0) r <- r + c(-0.5, 0.5)
  r
}

.fc_axis_ticks <- function(lim, n = 8L, special = NULL) {
  ticks <- pretty(lim, n = n)
  ticks <- ticks[ticks >= lim[1L] & ticks <= lim[2L]]
  special <- special[is.finite(special)]
  special <- special[special >= lim[1L] & special <= lim[2L]]
  ticks <- c(ticks, special)
  if (length(ticks) == 0L) ticks <- lim
  sort(unique(ticks))
}

.fc_axis_labels <- function(ticks) {
  ticks <- as.numeric(ticks)
  finite <- ticks[is.finite(ticks)]
  nonzero <- abs(finite[finite != 0])
  scientific <- length(finite) > 0L && (
    max(abs(finite)) >= 1e6 ||
      (length(nonzero) > 0L && min(nonzero) < 1e-3)
  )
  if (scientific) {
    return(format(ticks, trim = TRUE, scientific = TRUE, digits = 3))
  }
  vapply(ticks, function(x) {
    if (!is.finite(x)) return("")
    if (abs(x - round(x)) < 1e-8) {
      return(format(round(x), trim = TRUE, scientific = FALSE))
    }
    sub("\\.?0+$", "", format(round(x, 2L), trim = TRUE,
                              scientific = FALSE))
  }, character(1L))
}

.fc_text_width <- function(labels, cex = 1) {
  labels <- labels[nzchar(labels)]
  if (length(labels) == 0L) return(0)
  out <- tryCatch(
    graphics::strwidth(labels, units = "inches", cex = cex),
    error = function(e) rep(0.11 * cex, length(labels))
  )
  max(out, na.rm = TRUE)
}

.fc_text_height <- function(labels, cex = 1) {
  labels <- labels[nzchar(labels)]
  if (length(labels) == 0L) return(0)
  out <- tryCatch(
    graphics::strheight(labels, units = "inches", cex = cex),
    error = function(e) rep(0.16 * cex, length(labels))
  )
  max(out, na.rm = TRUE)
}

.fc_line_height <- function() {
  line <- graphics::par("cin")[2L] * graphics::par("lheight")
  if (!is.finite(line) || line <= 0) line <- 0.2
  line
}

.fc_trace_layout <- function(x_labels, y_labels, xlab, ylab, main) {
  cex.axis <- 0.95
  cex.lab <- 1.00
  cex.main <- 1.05
  line <- .fc_line_height()

  y_tick_width <- .fc_text_width(y_labels, cex.axis)
  y_title_width <- .fc_text_height(ylab, cex.lab)
  x_tick_height <- .fc_text_height(x_labels, cex.axis)
  x_title_height <- .fc_text_height(xlab, cex.lab)
  main_height <- .fc_text_height(main, cex.main)

  axis_label_line <- 0.75
  xlab_line <- axis_label_line + x_tick_height / line + 1.25
  ylab_line <- axis_label_line + y_tick_width / line + 1.25
  bottom <- xlab_line + x_title_height / line + 0.75
  left <- ylab_line + y_title_width / line + 0.85
  top <- main_height / line + 1.8
  right <- 1.2

  graphics::par(
    mar = c(bottom, left, top, right),
    mgp = c(2.0, axis_label_line, 0),
    tcl = -0.25
  )
  list(cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
       xlab.line = xlab_line, ylab.line = ylab_line)
}

.fc_protected_ticks <- function(ticks, protected) {
  if (length(protected) == 0L) return(rep(FALSE, length(ticks)))
  tol <- max(1e-8, diff(range(ticks, finite = TRUE)) * 1e-10)
  vapply(ticks, function(x) any(abs(x - protected) <= tol), logical(1L))
}

.fc_select_x_ticks <- function(ticks, labels, xlim, protected = NULL,
                               cex.axis = 1) {
  ord <- order(ticks)
  ticks <- ticks[ord]
  labels <- labels[ord]
  protected <- protected[is.finite(protected)]
  is.protected <- .fc_protected_ticks(ticks, protected)

  plot_width <- graphics::par("pin")[1L]
  if (!is.finite(plot_width) || plot_width <= 0) plot_width <- 5
  label_width <- tryCatch(
    graphics::strwidth(labels, units = "inches", cex = cex.axis),
    error = function(e) rep(0.10 * cex.axis, length(labels))
  )
  label_width <- label_width / plot_width * diff(xlim)
  pad <- 0.015 * diff(xlim)

  keep <- rep(FALSE, length(ticks))
  for (i in seq_along(ticks)) {
    overlap <- which(keep & abs(ticks[i] - ticks) <
                       (label_width[i] + label_width) / 2 + pad)
    if (length(overlap) == 0L) {
      keep[i] <- TRUE
    } else if (is.protected[i]) {
      drop <- overlap[!is.protected[overlap]]
      keep[drop] <- FALSE
      keep[i] <- TRUE
    }
  }

  list(ticks = ticks[keep], labels = labels[keep])
}

.fc_stan_loglik_name <- function(stan.obj) {
  mcmc_names <- tryCatch(names(stan.obj@sim$samples[[1L]]),
                         error = function(e) NULL)
  if (!is.null(mcmc_names)) {
    base_mcmc <- sub("\\[.*$", "", mcmc_names)
    if ("log_lik" %in% base_mcmc) return("log_lik")
    if ("log_lik_group" %in% base_mcmc) return("log_lik_group")
  }

  pars <- tryCatch(
    rstan::summary(stan.obj)$summary,
    error = function(e) NULL
  )
  if (is.null(pars)) {
    stop("Cannot read Stan summary from 'stan.obj'.", call. = FALSE)
  }
  base <- sub("\\[.*$", "", rownames(pars))
  if ("log_lik" %in% base) return("log_lik")
  if ("log_lik_group" %in% base) return("log_lik_group")
  stop("Stan object does not contain 'log_lik' or 'log_lik_group'.",
       call. = FALSE)
}

.fc_stan_loglik_weights <- function(x, par_name, n_loglik) {
  if (par_name != "log_lik_group" || !inherits(x, "FCGDINA")) {
    return(rep(1, n_loglik))
  }
  prep <- fcgdina_prepare(
    data = x$response,
    Q.matrix = x$Q.matrix,
    model = x$model,
    block.items = x$block.items,
    fc.type = x$fc.type
  )
  if (length(prep$response.count) != n_loglik) {
    return(rep(1, n_loglik))
  }
  prep$response.count
}

.fc_stan_loglik_trace <- function(x) {
  par_name <- .fc_stan_loglik_name(x$stan.obj)
  arr <- tryCatch(
    rstan::extract(
      x$stan.obj, pars = par_name, permuted = FALSE, inc_warmup = TRUE
    ),
    error = function(e) {
      stop("Cannot extract Stan log-likelihood draws: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  d <- dim(arr)
  if (is.null(d) || length(d) < 2L) {
    stop("Stan log-likelihood draws have an unsupported shape.",
         call. = FALSE)
  }
  if (length(d) == 2L) {
    out <- arr
  } else {
    n_loglik <- d[3L]
    weights <- .fc_stan_loglik_weights(x, par_name, n_loglik)
    out <- matrix(0, nrow = d[1L], ncol = d[2L])
    for (j in seq_len(n_loglik)) {
      slice <- arr[, , j, drop = FALSE]
      dim(slice) <- d[1:2]
      out <- out + slice * weights[j]
    }
  }
  storage.mode(out) <- "double"
  colnames(out) <- paste0("chain ", seq_len(ncol(out)))
  out
}

.fc_stan_stable_boundary <- function(stan.obj, n_iter) {
  sim <- stan.obj@sim
  iter <- sim$iter %||% n_iter
  warmup <- sim$warmup %||% 0L
  thin <- max(sim$thin %||% 1L)
  if (!is.numeric(warmup) || length(warmup) == 0L ||
      !is.finite(max(warmup)) || max(warmup) <= 0L ||
      !is.numeric(iter) || length(iter) == 0L || !is.finite(max(iter)) ||
      !is.finite(thin) || thin < 1L) {
    return(NA_real_)
  }
  post_warmup_saved <- if (!is.null(sim$n_save) && !is.null(sim$warmup2)) {
    max(as.numeric(sim$n_save) - as.numeric(sim$warmup2), na.rm = TRUE)
  } else {
    ceiling((max(iter) - max(warmup)) / thin)
  }
  boundary <- n_iter - post_warmup_saved + 1L
  if (!is.finite(boundary) || boundary <= 1 || boundary > n_iter) {
    return(NA_real_)
  }
  boundary
}

.fc_istem_loglik_trace <- function(x) {
  trace <- x$iStEM$logLik.trace
  if (is.null(trace)) {
    stop("iStEM object does not contain 'iStEM$logLik.trace'. Refit the ",
         "model to record iteration-level log-likelihood values.",
         call. = FALSE)
  }
  if (is.data.frame(trace)) {
    iter <- trace$iter %||% trace$batch %||% seq_len(nrow(trace))
    value <- trace$logLik %||% trace$log_lik %||% trace$loglik
  } else {
    value <- as.numeric(trace)
    iter <- seq_along(value)
  }
  value <- as.numeric(value)
  iter <- as.numeric(iter)
  keep <- is.finite(iter) & is.finite(value)
  if (!any(keep)) {
    stop("iStEM log-likelihood trace has no finite values.", call. = FALSE)
  }
  data.frame(iter = iter[keep], logLik = value[keep])
}

.fc_istem_boundary <- function(x) {
  boundary <- x$iStEM$stable.start
  if (is.null(boundary)) {
    boundary <- x$iStEM$stable.start.batch
  }
  if (!is.numeric(boundary) || length(boundary) != 1L ||
      !is.finite(boundary) || boundary <= 1) {
    return(NA_real_)
  }
  boundary
}

.fc_em_loglik_trace <- function(x) {
  trace <- x$EM$logLik.trace
  if (is.null(trace)) {
    stop("EM object does not contain 'EM$logLik.trace'.", call. = FALSE)
  }
  iter <- trace$iter %||% seq_len(nrow(trace))
  value <- trace$logLik %||% trace$log_lik %||% trace$loglik
  value <- as.numeric(value)
  iter <- as.numeric(iter)
  keep <- is.finite(iter) & is.finite(value)
  if (!any(keep)) {
    stop("EM log-likelihood trace has no finite values.", call. = FALSE)
  }
  data.frame(iter = iter[keep], logLik = value[keep])
}

.fc_plot_trace_frame <- function(xlim, ylim, xlab, ylab, main,
                                 boundary = NA_real_,
                                 legend_labels = NULL) {
  protected <- c(xlim[2L], boundary)
  x_ticks <- .fc_axis_ticks(xlim, n = 9L, special = protected)
  y_ticks <- .fc_axis_ticks(ylim)
  x_labels <- .fc_axis_labels(x_ticks)
  y_labels <- .fc_axis_labels(y_ticks)
  layout <- .fc_trace_layout(
    x_labels = x_labels, y_labels = y_labels, xlab = xlab, ylab = ylab,
    main = main
  )
  x_axis <- .fc_select_x_ticks(
    x_ticks, x_labels, xlim, protected = protected,
    cex.axis = layout$cex.axis
  )
  graphics::plot(
    NA_real_, NA_real_, type = "n", xlim = xlim, ylim = ylim,
    xlab = "", ylab = "", main = "",
    bty = "l", col.axis = "#374151", col.lab = "#374151",
    col.main = "#1F2937", axes = FALSE
  )
  graphics::grid(col = "#ECEFF1", lty = 1)
  graphics::axis(1, at = x_axis$ticks, labels = x_axis$labels,
                 col.axis = "#374151", cex.axis = layout$cex.axis)
  graphics::axis(2, at = y_ticks, labels = y_labels,
                 las = 1, col.axis = "#374151", cex.axis = layout$cex.axis)
  graphics::mtext(xlab, side = 1L, line = layout$xlab.line,
                  col = "#374151", cex = layout$cex.lab)
  graphics::mtext(ylab, side = 2L, line = layout$ylab.line,
                  col = "#374151", cex = layout$cex.lab)
  graphics::mtext(main, side = 3L, line = 0.6, font = 2L,
                  col = "#1F2937", cex = layout$cex.main)
  graphics::box(bty = "l")
  invisible(layout)
}

.fc_draw_boundary <- function(boundary, y_limits) {
  if (!is.finite(boundary)) return(invisible(NULL))
  graphics::abline(v = boundary, col = "#6B7280", lty = 2, lwd = 1.3)
  invisible(NULL)
}

.fc_trace_data <- function(x) {
  cls <- class(x)[1L]
  if (!is.null(x$stan.obj)) {
    logLik <- .fc_stan_loglik_trace(x)
    n_iter <- nrow(logLik)
    return(list(
      iter = seq_len(n_iter),
      logLik = logLik,
      boundary = .fc_stan_stable_boundary(x$stan.obj, n_iter),
      xlab = "Stan iteration",
      main = paste(cls, "Stan logLik trace"),
      legend_labels = if (ncol(logLik) > 1L) colnames(logLik) else NULL
    ))
  }

  if (!is.null(x$iStEM) && length(x$iStEM) > 0L) {
    trace <- .fc_istem_loglik_trace(x)
    return(list(
      iter = trace$iter,
      logLik = matrix(trace$logLik, ncol = 1L),
      boundary = .fc_istem_boundary(x),
      xlab = "iStEM iteration",
      main = paste(cls, "iStEM logLik trace"),
      legend_labels = NULL
    ))
  }

  if (!is.null(x$EM) && length(x$EM) > 0L) {
    trace <- .fc_em_loglik_trace(x)
    return(list(
      iter = trace$iter,
      logLik = matrix(trace$logLik, ncol = 1L),
      boundary = NA_real_,
      xlab = "EM iteration",
      main = paste(cls, "EM logLik trace"),
      legend_labels = NULL
    ))
  }

  stop("plot() requires Stan, iStEM, or EM log-likelihood trace data.",
       call. = FALSE)
}

.fc_draw_trace_legend <- function(labels, cols, boundary = NA_real_) {
  show_boundary <- is.finite(boundary)
  if (length(labels) == 0L && !show_boundary) return(invisible(NULL))
  legend_labels <- c(labels, if (show_boundary) "stable chain starts")
  legend_cols <- c(cols[seq_along(labels)], if (show_boundary) "#6B7280")
  legend_lty <- c(rep(1, length(labels)), if (show_boundary) 2)
  legend_lwd <- c(rep(2.0, length(labels)), if (show_boundary) 1.3)
  graphics::legend(
    "bottomright", legend = legend_labels, col = legend_cols,
    lty = legend_lty, lwd = legend_lwd, bty = "n",
    cex = max(0.7, min(0.85, 1.4 / sqrt(length(legend_labels)))),
    inset = 0.02, xpd = FALSE
  )
  invisible(NULL)
}

.fc_draw_loglik_trace <- function(trace) {
  logLik <- as.matrix(trace$logLik)
  if (nrow(logLik) != length(trace$iter)) {
    stop("Log-likelihood trace length is inconsistent with iterations.",
         call. = FALSE)
  }
  n_line <- ncol(logLik)
  ylim <- .fc_finite_ylim(logLik)
  xlim <- .fc_finite_xlim(trace$iter)
  legend_labels <- trace$legend_labels

  .fc_plot_trace_frame(
    xlim, ylim, xlab = trace$xlab,
    ylab = "Log-likelihood", main = trace$main,
    boundary = trace$boundary,
    legend_labels = legend_labels
  )
  .fc_draw_boundary(trace$boundary, ylim)

  cols <- .fc_trace_palette(n_line)
  lwd <- if (n_line > 1L) 2.0 else 1.8
  for (i in seq_len(n_line)) {
    graphics::lines(trace$iter, logLik[, i], col = cols[i], lwd = lwd)
  }
  .fc_draw_trace_legend(legend_labels, cols, boundary = trace$boundary)
  invisible(NULL)
}

.fc_plot_trace <- function(x, ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  .fc_draw_loglik_trace(.fc_trace_data(x))
}

.fc_plot_method <- function(x, y = NULL, ...) .fc_plot_trace(x, ...)

# ===========================================================================
# Plot methods
# ===========================================================================

#' @describeIn plot Log-likelihood trace plot for \code{MIRT} objects.
#' @method plot MIRT
#' @export
plot.MIRT <- .fc_plot_method

#' @describeIn plot Log-likelihood trace plot for \code{MGPCM} objects.
#' @method plot MGPCM
#' @export
plot.MGPCM <- .fc_plot_method

#' @describeIn plot Log-likelihood trace plot for \code{MGGUM} objects.
#' @method plot MGGUM
#' @export
plot.MGGUM <- .fc_plot_method

#' @describeIn plot Log-likelihood trace plot for \code{FCMIRT} objects.
#' @method plot FCMIRT
#' @export
plot.FCMIRT <- .fc_plot_method

#' @describeIn plot Log-likelihood trace plot for \code{FCDCM} objects.
#' @method plot FCDCM
#' @export
plot.FCDCM <- .fc_plot_method

#' @describeIn plot Log-likelihood trace plot for \code{FCGGUM} objects.
#' @method plot FCGGUM
#' @export
plot.FCGGUM <- .fc_plot_method

#' @describeIn plot Log-likelihood trace plot for \code{FCGDINA} objects.
#' @method plot FCGDINA
#' @export
plot.FCGDINA <- .fc_plot_method

#' @describeIn plot Log-likelihood trace plot for \code{TIRT} objects.
#' @method plot TIRT
#' @export
plot.TIRT <- .fc_plot_method
