## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----eval=FALSE---------------------------------------------------------------
# library(ForceChoice)
# 
# # Simulate binary response data
# sim <- sim.data.MIRT(N = 20, I = 6, D = 2, model = "2PL")
# 
# # Fit via iStEM
# fit <- fit.MIRT(sim$data, model = "2PL", D = 2,
#                 method = "iStEM",
#                 control.method = list(
#                   vis = FALSE, seed = 123,
#                   M = 2, B = 2, burnin.maxitr = 2,
#                   maxitr = 3, eps1 = 10, eps2 = 10,
#                   estimate.se = FALSE))
# 
# # Examine results
# print(fit)
# summary(fit)
# 
# # Item parameter estimates (first 6 items)
# head(coef(fit))
# 
# # Factor correlation matrix
# fit$Corr$est
# 
# # Trait recovery
# diag(cor(fit$theta$est, sim$theta))

## ----eval=FALSE---------------------------------------------------------------
# # Compute comprehensive fit indices
# gof <- get.fit.index(fit)
# 
# # Summary of fit indices
# summary(gof)
# 
# # Extract specific indices
# gof$M2       # Limited-information M2 statistic
# gof$RMSEA    # RMSEA with 90% CI
# gof$CFI      # Comparative Fit Index
# gof$TLI      # Tucker-Lewis Index
# gof$SRMSR    # Standardized Root Mean Square Residual
# gof$AIC      # Akaike Information Criterion
# gof$BIC      # Bayesian Information Criterion

## ----eval=FALSE---------------------------------------------------------------
# # Simulate forced-choice ranking data
# sim <- sim.data.FCMIRT(N.person = 20, N.block = 3, I.block = 2,
#                        D = 2, model = "2PL", fc.type = "RANK")
# 
# # The data contains ranking strings
# head(sim$data)
# 
# # Fit: block.items and fc.type are auto-detected
# fit <- fit.FCMIRT(sim$data, model = "2PL", D = 2,
#                   method = "iStEM",
#                   control.method = list(
#                     vis = FALSE, seed = 123,
#                     M = 2, B = 2, burnin.maxitr = 2,
#                     maxitr = 3, eps1 = 10, eps2 = 10,
#                     estimate.se = FALSE))
# 
# # Trait recovery
# cor(fit$theta$est, sim$theta)
# 
# # Goodness-of-fit (uses nominal binary expansion)
# gof <- get.fit.index(fit)
# summary(gof)

## ----eval=FALSE---------------------------------------------------------------
# sim <- sim.data.TIRT(N.person = 20, N.block = 3, I.block = 2,
#                      D = 2, fc.type = "RANK")
# 
# fit <- fit.TIRT(sim$data, Q.matrix = sim$Q.matrix,
#                 block.items = sim$block.items,
#                 method = "iStEM",
#                 control.method = list(
#                   vis = FALSE, seed = 123,
#                   M = 2, B = 2, burnin.maxitr = 2,
#                   maxitr = 3, eps1 = 10, eps2 = 10,
#                   estimate.se = FALSE))
# 
# # Structural parameters: loadings and uniquenesses
# head(coef(fit))
# 
# # Gamma matrix (pairwise intercepts)
# fit$gamma.matrix$est[1:5, 1:5]

## ----eval=FALSE---------------------------------------------------------------
# sim <- sim.data.FCDCM(N.person = 20, N.block = 3, D = 2,
#                       dcm.type = "DINA")
# 
# fit <- fit.FCDCM(sim$data, Q.matrix = sim$Q.matrix,
#                  block.items = sim$block.items,
#                  method = "iStEM",
#                  control.method = list(
#                    vis = FALSE, seed = 123,
#                    M = 2, B = 2, burnin.maxitr = 2,
#                    maxitr = 3, eps1 = 10, eps2 = 10,
#                    estimate.se = FALSE))
# 
# # Posterior attribute mastery probabilities
# head(fit$alpha$prob)
# 
# # Attribute mastery proportions
# colMeans(fit$alpha$prob > 0.5)
# 
# # Higher-order IRT parameters
# fit$delta$est

## ----eval=FALSE---------------------------------------------------------------
# sim <- sim.data.FCGDINA(N.person = 20, N.block = 2, I.block = 2,
#                         D = 2, model = "GDINA", fc.type = "RANK")
# 
# fit <- fit.FCGDINA(sim$data, Q.matrix = sim$Q.matrix,
#                    block.items = sim$block.items, model = "GDINA",
#                    fc.type = sim$fc.type, method = "EM",
#                    control.method = list(vis = FALSE, seed = 123,
#                                          maxitr = 2,
#                                          estimate.se = FALSE))
# 
# # Posterior attribute mastery probabilities
# head(fit$alpha$est)
# 
# # CDM item-parameter estimates
# coef(fit, type = "delta")

## ----eval=FALSE---------------------------------------------------------------
# sim <- sim.data.MGPCM(N = 20, I = 6, D = 2, length.poly = 4)
# 
# fit <- fit.MGPCM(sim$data, D = 2, method = "iStEM",
#                  control.method = list(
#                    vis = FALSE, seed = 123,
#                    M = 2, B = 2, burnin.maxitr = 2,
#                    maxitr = 3, eps1 = 10, eps2 = 10,
#                    estimate.se = FALSE))
# 
# # Category threshold parameters
# coef(fit)

## ----eval=FALSE---------------------------------------------------------------
# sim <- sim.data.MGGUM(N = 20, I = 6, D = 2, length.poly = 4)
# 
# fit <- fit.MGGUM(sim$data, D = 2, method = "iStEM",
#                  control.method = list(
#                    vis = FALSE, seed = 123,
#                    M = 2, B = 2, burnin.maxitr = 2,
#                    maxitr = 3, eps1 = 10, eps2 = 10,
#                    estimate.se = FALSE))
# 
# coef(fit)

## ----eval=FALSE---------------------------------------------------------------
# fit_rot <- rotate.MIRT(fit, rotate = "oblimin")
# 
# # Compare original and rotated loadings
# head(coef(fit))
# head(coef(fit_rot))

## ----eval=FALSE---------------------------------------------------------------
# # Item characteristic curves (ICC)
# plot(fit, type = "icc", items = 1:8)
# 
# # Person parameter distributions
# plot(fit, type = "theta")
# 
# # iStEM convergence trace
# plot(fit, type = "trace")

## ----eval=FALSE---------------------------------------------------------------
# fit <- fit.MIRT(data, model = "2PL", D = 2, method = "iStEM",
#                 control.method = list(
#                   vis = FALSE,
#                   seed = 123,      # Reproducibility
#                   M = 2, B = 2, burnin.maxitr = 2,
#                   maxitr = 3, eps1 = 10, eps2 = 10,
#                   estimate.se = FALSE))

## ----eval=FALSE---------------------------------------------------------------
# # stan code, long time
# fit <- fit.MIRT(data, model = "2PL", D = 2, method = "stan",
#                 control.method = list(chains = 1, iter = 200,
#                                       warmup = 100, cores = 1,
#                                       seed = 123))

## ----eval=FALSE---------------------------------------------------------------
# fit <- fit.MIRT(data, model = "2PL", D = 2, method = "iStEM",
#                 control.method = list(
#                   seed = 123,
#                   M = 2,            # Burn-in batches for convergence
#                   B = 2,            # Iterations per batch
#                   burnin.maxitr = 2,
#                   maxitr = 3,
#                   eps1 = 1.5,       # Geweke convergence threshold
#                   eps2 = 0.4        # MC error tolerance
#                 ))

