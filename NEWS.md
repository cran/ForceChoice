# ForceChoice 1.0.0 (2026-07-12)

## New Features

### Eight Model Families
- `fit.MIRT()` - Multidimensional Item Response Theory (1PL--4PL) for binary
  dominance data, with optional Q-matrix structure.
- `fit.MGPCM()` - Multidimensional Generalized Partial Credit Model for
  polytomous ordered-category responses.
- `fit.MGGUM()` - Multidimensional Generalized Graded Unfolding Model for
  ideal-point polytomous responses with signed Q-matrix.
- `fit.FCMIRT()` - Forced-Choice MIRT with item-level dominance endorsement
  and Luce--Plackett block-level ranking (RANK, MOLE, PICK).
- `fit.FCGGUM()` - Forced-Choice GGUM with item-level ideal-point endorsement
  and block-level ranking.
- `fit.TIRT()` - Thurstonian IRT for forced-choice with pairwise probit
  comparisons and latent utility differences.
- `fit.FCDCM()` - Forced-Choice Diagnostic Classification Model with
  higher-order latent trait, DINA/DINO condensation rules, and exact
  attribute-profile marginalization.
- `fit.FCGDINA()` - Forced-Choice GDINA model with DINA, DINO, ACDM, and
  GDINA item-level structures for ranking, most-least, and pick responses.

### Estimation Backends
- **Stan** (`method = "stan"`): Full Bayesian inference via Hamiltonian Monte
  Carlo (NUTS/HMC) with `rstan`.
- **iStEM** (`method = "iStEM"`): Fast iterative Stochastic EM with
  Metropolis-within-Gibbs person sampling and L-BFGS-B item optimization.
- **EM** (`method = "EM"`): Deterministic posterior-weight EM for FCGDINA.

### Comprehensive Goodness-of-Fit
- Limited-information M2 statistic (Maydeu-Olivares & Joe, 2005, 2006)
- RMSEA with confidence intervals, SRMSR, RMSR
- Comparative fit: CFI, TLI, IFI
- Information criteria: AIC, AICc, BIC, CAIC, SABIC, HQIC
- Pseudo-R2 measures: McFadden, Cox--Snell, Nagelkerke, and more
- Local dependence: Yen's Q3 (raw and adjusted)
- Posterior diagnostics: entropy, classification accuracy

### Data Simulation
- `sim.data.MIRT()`, `sim.data.MGPCM()`, `sim.data.MGGUM()` for traditional
  item response data.
- `sim.data.FCMIRT()`, `sim.data.FCGGUM()`, `sim.data.TIRT()`,
  `sim.data.FCDCM()`, `sim.data.FCGDINA()` for forced-choice data.
- All simulators support controlled parameter generation, factor correlation,
  and rotation.

### Solution Rotation
- `rotate.MIRT()` and `rotate.matrix()` for post-hoc rotation of MIRT
  solutions using promax or any GPArotation method.

### S3 Methods
- Consistent generics across all model classes: `coef()`,
  `confint()`, `deviance()`, `fitted()`, `logLik()`, `nobs()`, `plot()`,
  `predict()`, `print()`, `residuals()`, `summary()`, `update()`, `vcov()`.

## Bug Fixes

- This is the initial CRAN release.

## Notes

- The package was previously developed under the name **fcirt** and has been
  substantially rewritten and expanded to support eight model families and the
  iStEM estimation backend.
