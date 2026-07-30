# library(ForceChoice)
#
# ############################### data ###############################
# seed <- 13549
# N <- 500
# I <- 10
# D <- 2
# length.poly <- 3
# Q.matrix <- NULL
# Corr <- NULL
# rotate <- NULL
#
# set.seed(seed)
# data.obj <- sim.data.MGPCM(N = N, I = I, D = D, length.poly = length.poly,
#                            Q.matrix = Q.matrix, Corr = Corr, rotate = rotate)
#
# response <- data.obj$response
# head(response)
#
# # ############################### Stan MCMC ###############################
# # object.MGPCM.stan <- fit.MGPCM(response, D = D, Q.matrix = Q.matrix,
# #                                method = "stan",
# #                                control.method = c(list(
# #                                  init = "random",
# #                                  algorithm = "HMC"), list(seed = seed, vis = TRUE)))
# # cbind(data.obj$par, object.MGPCM.stan$par$est)
# # get.fit.index(object.MGPCM.stan)
#
# ############################### iStEM ###############################
# iStEM.control.method <- list(
#   burnin.maxitr = 40,
#   maxitr = 500,
#   theta.proposal.sd = 0.35,
#   optim.maxit = 50
# )
#
# object.MGPCM.iStEM <- fit.MGPCM(response, D = D, Q.matrix = Q.matrix,
#                                 method = "iStEM",
#                                 control.method = c(iStEM.control.method, list(seed = seed, vis = TRUE)))
# cbind(data.obj$par, object.MGPCM.iStEM$par$est)
#
# # ---- Validate logLik with new quadrature helpers ----
# ll <- logLik(object.MGPCM.iStEM)
# stopifnot(inherits(ll, "logLik"))
# stopifnot(is.finite(as.numeric(ll)))
# cat(sprintf("  MGPCM logLik: %.2f (df=%d)\n", as.numeric(ll), attr(ll, "df")))
#
# get.fit.index(object.MGPCM.iStEM)
