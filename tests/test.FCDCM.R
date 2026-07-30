# library(ForceChoice)
#
# ############################### data ###############################
# N.person <- 1000
# N.block <- 30
# D <- 4
# dcm.type <- "DINA"
#
# seed <- 23548
# set.seed(seed)
# data.obj <- sim.data.FCDCM(N.person = N.person, N.block = N.block,
#                            D = D, dcm.type = dcm.type)
# response <- data.obj$response
# data <- data.obj$data
# Q.matrix <- data.obj$Q.matrix
# block.items <- data.obj$block.items
#
# data.from.response <- get.data.from.response.FCDCM(response,
#                                                                  block.items)
# response.obj <- get.response.from.data.FCDCM(data, block.items)
# stopifnot(all(as.matrix(data) == as.matrix(data.from.response)))
# stopifnot(all(response == response.obj$response))
# block.items.detected <- get.block.items.from.data.FCDCM(data)
# stopifnot(all(vapply(seq_along(block.items), function(b) {
#   setequal(block.items[[b]], block.items.detected[[b]])
# }, logical(1L))))
#
# # ############################### Stan MCMC ###############################
# stan.control.method <- list(
#   init = "random",
#   algorithm = "HMC",
#   seed = seed,
#   vis = TRUE
# )
# object.FCDCM.stan <- fit.FCDCM(data, Q.matrix = Q.matrix,
#                                block.items = block.items,
#                                dcm.type = dcm.type,
#                                method = "stan",
#                                control.method = stan.control.method)
# cbind(data.obj$par, object.FCDCM.stan$par$est)
# cbind(data.obj$delta1, data.obj$delta0)
# object.FCDCM.stan$delta$est
# object.FCDCM.stan$par$Rhat
# get.fit.index(object.FCDCM.stan)
# mean(abs((object.FCDCM.stan$alpha$est - data.obj$alpha)))
#
# ############################## S3 ###############################
# class(object.FCDCM.stan)
# methods(class = class(object.FCDCM.stan))
#
# coef(object.FCDCM.stan, "par")
# confint(object.FCDCM.stan)
# deviance(object.FCDCM.stan)
# extract(object.FCDCM.stan, "par")
# fitted(object.FCDCM.stan)
# get.fit.index(object.FCDCM.stan)
# logLik(object.FCDCM.stan)
# nobs(object.FCDCM.stan)
# plot(object.FCDCM.stan)
# predict(object.FCDCM.stan)
# residuals(object.FCDCM.stan)
# temp <- summary(object.FCDCM.stan)
# update(object.FCDCM.stan, control.model=list(L=31))
# vcov(object.FCDCM.stan)
#
# ############################### iStEM ###############################
# iStEM.control.method <- list(
#   burnin.maxitr = 100,
#   maxitr = 500,
#   optim.maxit = 100,
#   seed = seed,
#   vis = TRUE
# )
#
# object.FCDCM.iStEM <- fit.FCDCM(data, Q.matrix = Q.matrix,
#                                 block.items = block.items,
#                                 dcm.type = dcm.type,
#                                 method = "iStEM",
#                                 control.method = iStEM.control.method)
# cbind(data.obj$par, object.FCDCM.iStEM$par$est)
# cbind(data.obj$delta1, data.obj$delta0)
# object.FCDCM.iStEM$delta$est
# get.fit.index(object.FCDCM.iStEM)
#
# sqrt(mean((data.obj$par - object.FCDCM.iStEM$par$est)**2))
# 1-mean(abs((object.FCDCM.iStEM$alpha$est - data.obj$alpha)))
#
# # ############################## S3 ###############################
# class(object.FCDCM.iStEM)
# methods(class = class(object.FCDCM.iStEM))
#
# coef(object.FCDCM.iStEM, "par")
# confint(object.FCDCM.iStEM)
# deviance(object.FCDCM.iStEM)
# extract(object.FCDCM.iStEM, "par")
# fitted(object.FCDCM.iStEM)
# get.fit.index(object.FCDCM.iStEM)
# logLik(object.FCDCM.iStEM)
# nobs(object.FCDCM.iStEM)
# plot(object.FCDCM.iStEM)
# predict(object.FCDCM.iStEM)
# residuals(object.FCDCM.iStEM)
# temp <- summary(object.FCDCM.iStEM)
# update(object.FCDCM.iStEM, control.model=list(L=31))
# vcov(object.FCDCM.iStEM)
