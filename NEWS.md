# spatialHVI 0.1.0 (2026-11-17)

- Add `fit_spatial_hvi()` with the manuscript's frozen MFVI updates, model
  priors, initialization, discrete rho grid, and mean-inverse residual summary.
- Add original joint variational intervals, discrete-grid LRVB covariance,
  logit-delta intervals, and separate simulation calibration. Both individual
  traits and bilateral means retain the original MFVI point estimates.
- Include the full-precision S17 calibration coefficients and a function for
  learning new coefficients from separate development simulations.
- Correct the Cholesky multiplication order in the legacy PCH summary.
  Trait-wise legacy point estimates are unchanged.
- Add tests for generalized eigenvalues, pair permutation invariance, response
  derivatives, bilateral covariance, variational draw quantiles, input alignment,
  and calibration. Retain the original fitting interfaces for compatibility.

# spatialHVI 0.0.0.9000

- Initial development package for halves-ordered spatial heritability analysis.
