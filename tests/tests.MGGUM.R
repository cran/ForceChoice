# library(ForceChoice)
#
# ############################### data ###############################
# seed <- 4592166
# N <- 500
# I <- 15
# D <- 2
# length.poly <- 3
# Corr <- NULL
# Q.matrix <- NULL
#
# set.seed(seed)
# data.obj <- sim.data.MGGUM(N = N, I = I, D = D, length.poly = length.poly,
#                            Q.matrix = Q.matrix, Corr = Corr)
#
# data.obj$par
# response <- data.obj$response
# head(response)
# apply(response, 2, table)
# Q.matrix <- data.obj$arguments$Q.matrix
#
# # ############################### Stan MCMC ###############################
# # object.MGGUM.stan <- fit.MGGUM(response, D = D, Q.matrix = Q.matrix,
# #                                length.poly = length.poly,
# #                                method = "stan",
# #                                control.method = c(list(
# #                                  init = "random",
# #                                  algorithm = "HMC"), list(seed = seed, vis = TRUE)))
# # cbind(data.obj$par, object.MGGUM.stan$par$est)
# # object.MGGUM.stan$par$Rhat
# # get.fit.index(object.MGGUM.stan)
#
# ############################### iStEM ###############################
# iStEM.control.method <- list(
#   burnin.maxitr = 40,
#   maxitr = 500,
#   theta.proposal.sd = 0.35,
#   optim.maxit = 50
# )
#
# object.MGGUM.iStEM <- fit.MGGUM(response, D = D, Q.matrix = Q.matrix,
#                                 length.poly = length.poly,
#                                 method = "iStEM",
#                                 control.method = c(iStEM.control.method, list(seed = seed, vis = TRUE)))
# cbind(data.obj$par, object.MGGUM.iStEM$par$est)
#
# # ---- Validate logLik with new quadrature helpers ----
# ll <- logLik(object.MGGUM.iStEM)
# stopifnot(inherits(ll, "logLik"))
# stopifnot(is.finite(as.numeric(ll)))
# cat(sprintf("  MGGUM logLik: %.2f (df=%d)\n", as.numeric(ll), attr(ll, "df")))
#
# get.fit.index(object.MGGUM.iStEM)
