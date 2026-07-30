# ForceChoice

## Forced-Choice Modeling Based on Item Response Theory and Cognitive Diagnosis

The **ForceChoice** package provides a unified framework for fitting,
simulating, and evaluating forced-choice and traditional item response
models. It supports **eight model families** spanning dominance,
ideal-point/unfolding, and diagnostic classification paradigms, with both
full Bayesian estimation (Stan/HMC) and a fast iterative stochastic EM
(iStEM) algorithm. FCGDINA also provides a deterministic EM estimator.

The model descriptions below distinguish the package's implemented
parameterizations from the broader source literature. Core references are
given with DOI links where available so that the modeling assumptions,
response-process assumptions, and goodness-of-fit framework can be traced
back to peer-reviewed sources.

## Features

- **Eight model families**: MIRT, MGPCM, MGGUM, FCMIRT, FCGGUM, TIRT,
  FCDCM, FCGDINA
- **Estimation backends**: Stan (full Bayesian HMC), iStEM (fast stochastic
  EM for large data), and deterministic EM for FCGDINA
- **Comprehensive goodness-of-fit**: Limited-information M2, RMSEA, CFI, TLI,
  SRMSR, pseudo-*R*^2^, Yen's Q3, posterior diagnostics
- **Forced-choice response types**: Full ranking (RANK), most-least (MOLE),
  and best-only (PICK)
- **Consistent S3 interface**: `coef()`, `fitted()`, `logLik()`,
  `plot()`, `predict()`, `print()`, `residuals()`, `summary()`, `update()`,
  `vcov()`, `confint()`

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Naidantu/ForceChoice")
```

## Quick Start

### Example 1: Multidimensional IRT (MIRT)

```r
library(ForceChoice)

# Simulate binary response data from a 2-dimensional 2PL model
sim <- sim.data.MIRT(N = 20, I = 6, D = 2, model = "2PL")

# Fit via iStEM (fast, scalable to large datasets)
fit <- fit.MIRT(sim$data, model = "2PL", D = 2, method = "iStEM",
                control.method = list(
                  vis = FALSE, seed = 123,
                  M = 2, B = 2, burnin.maxitr = 2,
                  maxitr = 3, eps1 = 10, eps2 = 10,
                  estimate.se = FALSE))

# Examine results
print(fit)
summary(fit)

# Item parameter estimates
head(coef(fit))

# Person trait estimates
head(fit$theta$est)

# Factor correlation
fit$Corr$est

# Goodness-of-fit evaluation
gof <- get.fit.index(fit)
summary(gof)
```

### Example 2: Forced-Choice MIRT (FCMIRT)

```r
# Simulate forced-choice ranking data
sim <- sim.data.FCMIRT(N.person = 20, N.block = 3, I.block = 2,
                       D = 2, model = "2PL", fc.type = "RANK")

# The block.items and fc.type are auto-detected from the data
fit <- fit.FCMIRT(sim$data, model = "2PL", D = 2,
                  method = "iStEM",
                  control.method = list(
                    vis = FALSE, seed = 123,
                    M = 2, B = 2, burnin.maxitr = 2,
                    maxitr = 3, eps1 = 10, eps2 = 10,
                    estimate.se = FALSE))

# Trait recovery: correlation between estimates and true values
diag(cor(fit$theta$est, sim$theta))

# Goodness-of-fit for forced-choice data
gof <- get.fit.index(fit)
summary(gof)
```

### Example 3: Thurstonian IRT (TIRT)

```r
# Simulate forced-choice data for TIRT
sim <- sim.data.TIRT(N.person = 20, N.block = 3, I.block = 2,
                     D = 2, fc.type = "RANK")

fit <- fit.TIRT(sim$data, Q.matrix = sim$Q.matrix,
                block.items = sim$block.items,
                method = "iStEM",
                control.method = list(
                  vis = FALSE, seed = 123,
                  M = 2, B = 2, burnin.maxitr = 2,
                  maxitr = 3, eps1 = 10, eps2 = 10,
                  estimate.se = FALSE))

# Trait recovery
cor(fit$theta$est, sim$theta)
```

### Example 4: Forced-Choice DCM (FCDCM)

```r
# Simulate paired-comparison data with DINA condensation rule
sim <- sim.data.FCDCM(N.person = 20, N.block = 3, D = 2,
                      dcm.type = "DINA")

# Fit
fit <- fit.FCDCM(sim$data, Q.matrix = sim$Q.matrix,
                 block.items = sim$block.items,
                 method = "iStEM",
                 control.method = list(
                   vis = FALSE, seed = 123,
                   M = 2, B = 2, burnin.maxitr = 2,
                   maxitr = 3, eps1 = 10, eps2 = 10,
                   estimate.se = FALSE))

# Posterior attribute mastery profiles
head(fit$alpha$est)

# Attribute mastery proportions
colMeans(fit$alpha$est > 0.5)
```

### Example 5: Forced-Choice GDINA (FCGDINA)

```r
# Simulate forced-choice diagnostic data under a CDM item model
sim <- sim.data.FCGDINA(N.person = 20, N.block = 2, I.block = 2,
                        D = 2, model = "GDINA", fc.type = "RANK")

# Fit with deterministic EM for a reproducible baseline
fit <- fit.FCGDINA(sim$data, Q.matrix = sim$Q.matrix,
                   block.items = sim$block.items, model = "GDINA",
                   fc.type = sim$fc.type, method = "EM",
                   control.method = list(vis = FALSE, seed = 123,
                                         maxitr = 2,
                                         estimate.se = FALSE))

head(fit$alpha$est)
```

## Model Overview

| Function | Model | Response Type | Key Reference |
|---|---|---|---|
| `fit.MIRT()` | Multidimensional IRT (1PL--4PL) | Binary | Reckase (2009) |
| `fit.MGPCM()` | Generalized Partial Credit | Polytomous | Muraki (1992) |
| `fit.MGGUM()` | Generalized Graded Unfolding | Polytomous | Roberts et al. (2000) |
| `fit.FCMIRT()` | Forced-Choice MIRT | Ranking/MOLE/PICK | Zheng et al. (2024); Luce (1959); Plackett (1975) |
| `fit.FCGGUM()` | Forced-Choice GGUM | Ranking/MOLE/PICK | Lee et al. (2018); Roberts et al. (2000) |
| `fit.TIRT()` | Thurstonian IRT | Ranking/MOLE/PICK | Brown & Maydeu-Olivares (2011) |
| `fit.FCDCM()` | Forced-Choice DCM | Paired comparison | Huang (2022) |
| `fit.FCGDINA()` | Forced-Choice GDINA | Ranking/MOLE/PICK | de la Torre (2011); Luce (1959); Plackett (1975) |

## Estimation Methods

### Stan (method = "stan")

Full Bayesian inference via Hamiltonian Monte Carlo (NUTS). Provides posterior
means, standard deviations, and R-hat convergence diagnostics. Recommended for
final inference and small-to-moderate datasets.

```r
# stan code, long time
fit <- fit.MIRT(data, model = "2PL", D = 2, method = "stan",
                control.method = list(chains = 1, iter = 200,
                                      warmup = 100, cores = 1,
                                      seed = 123))
```

### iStEM (method = "iStEM")

Iterative Stochastic EM with Metropolis-within-Gibbs person sampling and
L-BFGS-B item optimization. Scales to large datasets. Convergence monitored
via Geweke diagnostics and batch-means Monte Carlo error.

```r
fit <- fit.MIRT(data, model = "2PL", D = 2, method = "iStEM",
                control.method = list(
                  vis = FALSE, seed = 123,
                  M = 2, B = 2, burnin.maxitr = 2,
                  maxitr = 3, eps1 = 10, eps2 = 10,
                  estimate.se = FALSE))
```

### EM (method = "EM")

Deterministic posterior-weight EM is available for FCGDINA. It provides a
fast, reproducible point-estimation baseline without MCMC sampling.

## Goodness-of-Fit

All fitted models support comprehensive fit evaluation via `get.fit.index()`:

```r
gof <- get.fit.index(fit)
summary(gof)
```

Fit indices reported:
- **Information criteria**: AIC, AICc, BIC, CAIC, SABIC, HQIC
- **Absolute fit**: M2, RMSEA (with 90% CI), SRMSR
- **Comparative fit**: CFI, TLI, IFI
- **Pseudo-R2**: McFadden, Cox--Snell, Nagelkerke, and others
- **Local dependence**: Yen's Q3
- **Classification**: Posterior entropy, mean max posterior

## References

Brown, A., & Maydeu-Olivares, A. (2011). Item response modeling of
forced-choice questionnaires. *Educational and Psychological Measurement*,
*71*(3), 460--502. <https://doi.org/10.1177/0013164410375112>

de la Torre, J. (2011). The generalized DINA model framework.
*Psychometrika*, *76*(2), 179--199. <https://doi.org/10.1007/s11336-011-9207-7>

Huang, H.-Y. (2022). Diagnostic classification model for forced-choice
items and noncognitive tests. *Educational and Psychological Measurement*,
*83*(1), 146--180. <https://doi.org/10.1177/00131644211069906>

Lee, P., Joo, S.-H., Stark, S., & Chernyshenko, O. S. (2018). GGUM-RANK
statement and person parameter estimation with multidimensional forced
choice triplets. *Applied Psychological Measurement*, *43*(3), 226--240.
<https://doi.org/10.1177/0146621618768294>

Luce, R. D. (1959). *Individual choice behavior: A theoretical analysis*.
Wiley.

Maydeu-Olivares, A., & Joe, H. (2005). Limited- and full-information
estimation and goodness-of-fit testing in 2^n contingency tables: A unified
framework. *Journal of the American Statistical Association*, *100*(471),
1009--1020. <https://doi.org/10.1198/016214504000002069>

Maydeu-Olivares, A., & Joe, H. (2006). Limited information goodness-of-fit
testing in multidimensional contingency tables. *Psychometrika*, *71*(4),
713--732. <https://doi.org/10.1007/s11336-005-1295-9>

Muraki, E. (1992). A generalized partial credit model: Application of an EM
algorithm. *Applied Psychological Measurement*, *16*(2), 159--176.
<https://doi.org/10.1177/014662169201600206>

Plackett, R. L. (1975). The analysis of permutations. *Journal of the Royal
Statistical Society: Series C (Applied Statistics)*, *24*(2), 193--202.
<https://doi.org/10.2307/2346567>

Reckase, M. D. (2009). *Multidimensional Item Response Theory*. Springer.
<https://doi.org/10.1007/978-0-387-89976-3>

Roberts, J. S., Donoghue, J. R., & Laughlin, J. E. (2000). A general item
response theory model for unfolding unidimensional polytomous responses.
*Applied Psychological Measurement*, *24*(1), 3--32.
<https://doi.org/10.1177/01466216000241001>

Tu, N., Zhang, B., Angrave, L., Sun, T., & Neuman, M. (2023). Estimating the
multidimensional generalized graded unfolding model with covariates using a
Bayesian approach. *Journal of Intelligence*, *11*(8), 163.
<https://doi.org/10.3390/jintelligence11080163>

Zheng, C., Liu, J., Li, Y., Xu, P., Zhang, B., Wei, R., Zhang, W.,
Liu, B., & Huang, J. (2024). A 2PLM-RANK multidimensional forced-choice
model and its fast estimation algorithm. *Behavior Research Methods*,
*56*(6), 6363--6388. <https://doi.org/10.3758/s13428-023-02315-x>

Zhu, Y.-A., Xu, J., Wang, D., Li, X., Cai, Y., & Tu, D. (2024). A ranking
forced choice diagnostic classification model for psychological assessment
using forced choice questionnaires. *British Journal of Mathematical and
Statistical Psychology*, *78*(2), 617--646.
<https://doi.org/10.1111/bmsp.12376>

## License

GPL (>= 3)
