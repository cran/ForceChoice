# library(ForceChoice)
#
# ############################### data ###############################
# N.person <- 500L
# N.block <- 20
# I.block <- 2L
# D <- 2L
# model <- "GDINA"
# fc.type <- "RANK"
# control <- list(single=FALSE)
#
# seed <- 3546
# set.seed(seed)
#
# data.obj <- sim.data.FCGDINA(
#   N.person = N.person, N.block = N.block, I.block = I.block,
#   D = D, model = model, fc.type = fc.type, control = control
# )
#
# response <- data.obj$response
# data <- data.obj$data
# block.items <- data.obj$block.items
# Q.matrix <- data.obj$Q.matrix
# fc.type <- data.obj$fc.type
#
# # ---- Validate data conversion round-trip ----
# data.from.response <- get.data.from.response(
#   response, block.items, fc.type = fc.type
# )
# response.from.data <- get.response.from.data(
#   data, block.items, fc.type = fc.type
# )
# stopifnot(all(as.matrix(data) == as.matrix(data.from.response)))
# stopifnot(all(response == response.from.data))
#
# # ---- Verify unified block-item extraction ----
# block.items.from.data <- get.block.items.from.data(data)
# stopifnot(length(block.items.from.data) == N.block)
# stopifnot(all(mapply(setequal, block.items, block.items.from.data)))
#
# # ---- Verify pi estimation ----
# cat(sprintf("  True alpha patterns: %d\n", nrow(data.obj$alpha.patterns)))
# cat(sprintf("  Structural parameters (pi): %d free\n", 2^D - 1L))
#
# ############################## Stan MCMC ###############################
# stan.control.method <- list(
#   init   = "random",
#   algorithm = "HMC"
# )
#
# object.FCGDINA.stan <- fit.FCGDINA(data, Q.matrix = Q.matrix,
#                                     model = model,
#                                     block.items = block.items,
#                                     fc.type = fc.type,
#                                     method = "stan",
#                                     control.method = c(stan.control.method, list(seed = seed, vis = TRUE)))
# cbind(do.call(c, data.obj$delta), do.call(c, object.FCGDINA.stan$delta$est))
# object.FCGDINA.stan$delta$Rhat
# object.FCGDINA.stan$pi
# get.fit.index(object.FCGDINA.stan)
# mean(abs((object.FCGDINA.stan$alpha$est > 0.5) * 1L - data.obj$alpha))
#
# # ############################## S3 ###############################
# class(object.FCGDINA.stan)
# methods(class = class(object.FCGDINA.stan))
#
# coef(object.FCGDINA.stan, "par")
# confint(object.FCGDINA.stan)
# deviance(object.FCGDINA.stan)
# extract(object.FCGDINA.stan, "delta")
# fitted(object.FCGDINA.stan)
# get.fit.index(object.FCGDINA.stan)
# logLik(object.FCGDINA.stan)
# nobs(object.FCGDINA.stan)
# plot(object.FCGDINA.stan)
# predict(object.FCGDINA.stan)
# residuals(object.FCGDINA.stan)
# temp <- summary(object.FCGDINA.stan)
# update(object.FCGDINA.stan, control.model=list(L=31))
# vcov(object.FCGDINA.stan)
#
# ############################### EM ###############################
# EM.control.method <- list(
#   maxitr = 2000L,
#   minitr = 2L,
#   tol = 1e-4,
#   par.tol = 1e-4,
#   optim.maxit = 200L,
#   estimate.se = FALSE
# )
#
# object.FCGDINA.EM <- fit.FCGDINA(
#   data, Q.matrix = Q.matrix,
#   model = model,
#   block.items = block.items,
#   fc.type = fc.type,
#   control.method = c(EM.control.method, list(seed = seed, vis = TRUE))
# )
#
# cbind(do.call(c, data.obj$delta), do.call(c, object.FCGDINA.EM$delta$est))
# object.FCGDINA.EM$pi
# get.fit.index(object.FCGDINA.EM)
# mean(abs((object.FCGDINA.EM$alpha$est > 0.5) * 1L - data.obj$alpha))
#
# # ############################## S3 ###############################
# class(object.FCGDINA.EM)
# methods(class = class(object.FCGDINA.EM))
#
# coef(object.FCGDINA.EM, "par")
# confint(object.FCGDINA.EM)
# deviance(object.FCGDINA.EM)
# extract(object.FCGDINA.EM, "delta")
# fitted(object.FCGDINA.EM)
# get.fit.index(object.FCGDINA.EM)
# logLik(object.FCGDINA.EM)
# nobs(object.FCGDINA.EM)
# plot(object.FCGDINA.EM)
# predict(object.FCGDINA.EM)
# residuals(object.FCGDINA.EM)
# temp <- summary(object.FCGDINA.EM)
# # update(object.FCGDINA.EM, control.model=list(L=31))
# vcov(object.FCGDINA.EM)
#
# ############################### iStEM ###############################
# iStEM.control.method <- list(
#   burnin.maxitr = 40,
#   maxitr = 500
# )
#
# object.FCGDINA.iStEM <- fit.FCGDINA(
#   data, Q.matrix = Q.matrix,
#   model = model,
#   block.items = block.items,
#   fc.type = fc.type,
#   method = "iStEM",
#   control.method = c(iStEM.control.method, list(seed = seed, vis = TRUE))
# )
#
# cbind(do.call(c, data.obj$delta), do.call(c, object.FCGDINA.iStEM$delta$est))
# object.FCGDINA.iStEM$pi
# get.fit.index(object.FCGDINA.iStEM)
# mean(abs((object.FCGDINA.iStEM$alpha$est > 0.5) * 1L - data.obj$alpha))
#
# # ############################## S3 ###############################
# class(object.FCGDINA.iStEM)
# methods(class = class(object.FCGDINA.iStEM))
#
# coef(object.FCGDINA.iStEM, "par")
# confint(object.FCGDINA.iStEM)
# deviance(object.FCGDINA.iStEM)
# extract(object.FCGDINA.iStEM, "delta")
# fitted(object.FCGDINA.iStEM)
# get.fit.index(object.FCGDINA.iStEM)
# logLik(object.FCGDINA.iStEM)
# nobs(object.FCGDINA.iStEM)
# plot(object.FCGDINA.iStEM)
# predict(object.FCGDINA.iStEM)
# residuals(object.FCGDINA.iStEM)
# temp <- summary(object.FCGDINA.iStEM)
# # update(object.FCGDINA.iStEM, control.model=list(L=31))
# vcov(object.FCGDINA.iStEM)
