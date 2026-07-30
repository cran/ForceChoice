# library(ForceChoice)
#
# ############################### data ###############################
# N.person <- 500
# N.block <- 10
# I.block <- 3
# D <- 3
# fc.type <- "RANK"
# control <- NULL
#
# seed <- 3548875
# set.seed(seed)
# data.obj <- sim.data.TIRT(N.person = N.person, N.block = N.block,
#                           I.block = I.block, D = D, fc.type = fc.type,
#                           control = control)
# data.obj$par
# Q.matrix <- data.obj$Q.matrix
# response <- data.obj$response
# head(response)
# block.items <- data.obj$block.items
# pairs.value <- data.obj$pairs.value
# data <- get.data.from.response.TIRT(response, block.items,
#                                                   fc.type = fc.type,
#                                                   pairs.value = pairs.value)
# response.obj <- get.response.from.data.TIRT(data, block.items,
#                                                           fc.type = fc.type)
# head(response - response.obj$response)
# stopifnot(all(data == get.data.from.response.TIRT(response.obj$response,
#                                                       block.items,
#                                                       fc.type = fc.type,
#                                                       pairs.value = response.obj$pairs.value)))
#
# # ---- Verify unified block-item extraction ----
# block.items.from.data <- get.block.items.from.data(data)
# stopifnot(length(block.items.from.data) == N.block)
# stopifnot(all(mapply(setequal, block.items, block.items.from.data)))
#
# # ############################### Stan MCMC ###############################
# # object.TIRT.stan <- fit.TIRT(data, Q.matrix, block.items = block.items,
# #                              fc.type = fc.type,
# #                              method = "stan",
# #                              control.method = c(list(
# #                                init = "random",
# #                                algorithm = "HMC"), list(seed = seed, vis = TRUE)))
# # object.TIRT.stan
# # object.TIRT.stan$npar
# # cbind(data.obj$par, object.TIRT.stan$par$est)
# # object.TIRT.stan$par$Rhat
# # cbind(data.obj$gamma, object.TIRT.stan$gamma$est)
# # object.TIRT.stan$gamma$Rhat
# # get.fit.index(object.TIRT.stan)
#
# ############################### iStEM ###############################
# iStEM.control.method <- list(
#   burnin.maxitr = 40,
#   maxitr = 500,
#   theta.proposal.sd = 0.35,
#   optim.maxit = 50
# )
#
# object.TIRT.iStEM <- fit.TIRT(data, Q.matrix, block.items = block.items,
#                               fc.type = fc.type,
#                               method = "iStEM",
#                               control.method = c(iStEM.control.method, list(seed = seed, vis = TRUE)))
# object.TIRT.iStEM
# object.TIRT.iStEM$npar
# cbind(data.obj$par, object.TIRT.iStEM$par$est)
# cbind(data.obj$gamma, object.TIRT.iStEM$gamma$est)
#
# # ---- Validate logLik with new quadrature helpers ----
# ll <- logLik(object.TIRT.iStEM)
# stopifnot(inherits(ll, "logLik"))
# stopifnot(is.finite(as.numeric(ll)))
# cat(sprintf("  TIRT logLik: %.2f (df=%d)\n", as.numeric(ll), attr(ll, "df")))
#
# get.fit.index(object.TIRT.iStEM)
