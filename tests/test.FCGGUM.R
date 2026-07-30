# library(ForceChoice)
#
# ############################### data ###############################
# N.person <- 300
# N.block <- 10
# I.block <- 3
# D <- 3
# fc.type <- "MOLE"
# control <- NULL
#
# seed <- 459821
#
# set.seed(seed)
# data.obj <- sim.data.FCGGUM(N.person = N.person, N.block = N.block,
#                             I.block = I.block, D = D, fc.type = fc.type,
#                             control = control)
# data.obj$par
# response <- data.obj$response
# block.items <- data.obj$block.items
# Q.matrix <- data.obj$Q.matrix
# fc.type <- data.obj$fc.type
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
# # ---- Verify normalize_fc_type ----
# stopifnot(identical(
#   normalize_fc_type(NULL, 3L),
#   c("RANK", "RANK", "RANK")))
# stopifnot(identical(
#   normalize_fc_type(c("mole", "rank", "pick"), 3L),
#   c("MOLE", "RANK", "PICK")))
#
# # ############################### Stan MCMC ###############################
# # object.FCGGUM.stan <- fit.FCGGUM(data, Q.matrix = Q.matrix,
# #                                  block.items = block.items, D = D,
# #                                  fc.type = fc.type,
# #                                  method = "stan",
# #                                  control.method = c(list(
# #                                    init = "random",
# #                                    algorithm = "HMC"), list(seed = seed, vis = TRUE)))
# # cbind(data.obj$par, object.FCGGUM.stan$par$est)
# # object.FCGGUM.stan$par$Rhat
# # get.fit.index(object.FCGGUM.stan)
#
# ############################### iStEM ###############################
# iStEM.control.method <- list(
#   burnin.maxitr = 40,
#   maxitr = 500,
#   theta.proposal.sd = 0.35,
#   optim.maxit = 50
# )
#
# object.FCGGUM.iStEM <- fit.FCGGUM(data, Q.matrix = Q.matrix,
#                                   block.items = block.items, D = D,
#                                   fc.type = fc.type,
#                                   method = "iStEM",
#                                   control.method = c(iStEM.control.method, list(seed = seed, vis = TRUE)))
# cbind(data.obj$par, object.FCGGUM.iStEM$par$est)
#
# # ---- Validate logLik with new quadrature helpers ----
# ll <- logLik(object.FCGGUM.iStEM)
# stopifnot(inherits(ll, "logLik"))
# stopifnot(is.finite(as.numeric(ll)))
# stopifnot(!is.null(attr(ll, "theta.norm")))
# stopifnot(!is.null(attr(ll, "prob")))
# stopifnot(!is.null(attr(ll, "pi")))
# cat(sprintf("  FCGGUM logLik: %.2f (df=%d)\n", as.numeric(ll), attr(ll, "df")))
#
# get.fit.index(object.FCGGUM.iStEM)
