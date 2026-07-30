# library(ForceChoice)
#
# ############################### data ###############################
# N.person <- 500
# N.block <- 10
# I.block <- 2
# D <- 2
# model <- "m1pl"
# fc.type <- "MOLE"
# control <- NULL
#
# seed <- 362154
#
# set.seed(seed)
# data.obj <- sim.data.FCMIRT(N.person = N.person, N.block = N.block,
#                             I.block = I.block, D = D, model = model,
#                             fc.type = fc.type, control = control)
# data.obj$par
# response <- data.obj$response
# block.items <- data.obj$block.items
# Q.matrix <- data.obj$Q.matrix
# fc.type <- data.obj$fc.type
# block.b.sum <- vapply(block.items, function(items) {
#   sum(data.obj$par[items, D + 1L])
# }, numeric(1L))
# stopifnot(max(abs(block.b.sum)) < 1e-10)
# data <- get.data.from.response(response, block.items, fc.type = fc.type)
# head(data)
# head(response)
# all(data == data.obj$data)
# response.obj <- get.response.from.data(data, block.items, fc.type = fc.type)
# stopifnot(all(response == response.obj))
#
# # ---- Verify unified block-item extraction ----
# block.items.from.data <- get.block.items.from.data(data)
# stopifnot(length(block.items.from.data) == N.block)
# stopifnot(all(mapply(setequal, block.items, block.items.from.data)))
#
# # ############################### Stan MCMC ###############################
# # object.FCMIRT.stan <- fit.FCMIRT(data, model = model, Q.matrix = Q.matrix,
# #                                  block.items = block.items, D = D,
# #                                  fc.type = fc.type,
# #                                  method = "stan",
# #                                  control.method = c(list(
# #                                    init = "random",
# #                                    algorithm = "HMC"), list(seed = seed, vis = TRUE)))
# # cbind(data.obj$par, object.FCMIRT.stan$par$est)
# # object.FCMIRT.stan$par$Rhat
# # get.fit.index(object.FCMIRT.stan)
#
# ############################### iStEM ###############################
# iStEM.control.method <- list(
#   burnin.maxitr = 40,
#   maxitr = 500,
#   theta.proposal.sd = 0.10,
#   optim.maxit = 50
# )
#
# object.FCMIRT.iStEM <- fit.FCMIRT(data, model = model, Q.matrix = Q.matrix,
#                                   block.items = block.items, D = D,
#                                   fc.type = fc.type,
#                                   method = "iStEM",
#                                   control.method = c(iStEM.control.method, list(seed = seed, vis = TRUE)))
# cbind(data.obj$par, object.FCMIRT.iStEM$par$est)
#
# # ---- Validate logLik with new quadrature helpers ----
# ll <- logLik(object.FCMIRT.iStEM)
# stopifnot(inherits(ll, "logLik"))
# stopifnot(is.finite(as.numeric(ll)))
# stopifnot(!is.null(attr(ll, "theta.norm")))
# cat(sprintf("  FCMIRT logLik: %.2f (df=%d)\n", as.numeric(ll), attr(ll, "df")))
#
# get.fit.index(object.FCMIRT.iStEM)
