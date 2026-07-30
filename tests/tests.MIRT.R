# library(ForceChoice)
#
# ############################### data ###############################
# seed <- 4688567
# N <- 500
# I <- 10
# D <- 2
# Q.matrix <- NULL
# model <- "m2pl"
# Corr <- NULL
# rotate <- NULL
#
# set.seed(seed)
# data.obj <- sim.data.MIRT(N = N, I = I, D = D, model = model,
#                           Q.matrix = Q.matrix, Corr = Corr, rotate = rotate)
#
# response <- data.obj$response
# head(response)
#
# # ############################### Stan MCMC ###############################
# # object.MIRT.stan <- fit.MIRT(response, model = model, D = D,
# #                              Q.matrix = Q.matrix,
# #                              method = "stan",
# #                              control.method = c(list(
# #                                init = "random",
# #                                algorithm = "HMC"), list(seed = seed, vis = TRUE)))
# # cbind(data.obj$par, object.MIRT.stan$par$est)
# # get.fit.index(object.MIRT.stan)
#
# ############################### iStEM ###############################
# iStEM.control.method <- list(
#   burnin.maxitr = 40,
#   maxitr = 500,
#   theta.proposal.sd = 0.35,
#   optim.maxit = 50
# )
#
# object.MIRT.iStEM <- fit.MIRT(response, model = model, D = D,
#                               Q.matrix = Q.matrix,
#                               method = "iStEM",
#                               control.method = c(iStEM.control.method, list(seed = seed, vis = TRUE)))
# cbind(data.obj$par, object.MIRT.iStEM$par$est)
#
# # ---- Validate logLik with new quadrature helpers ----
# ll <- logLik(object.MIRT.iStEM)
# stopifnot(inherits(ll, "logLik"))
# stopifnot(is.finite(as.numeric(ll)))
# stopifnot(!is.null(attr(ll, "theta.norm")))
# cat(sprintf("  MIRT logLik: %.2f (df=%d)\n", as.numeric(ll), attr(ll, "df")))
#
# get.fit.index(object.MIRT.iStEM)
