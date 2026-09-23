# birmrssim 0.1.1

## Behavior changes

* `one_simulation()` now has an explicit argument list instead of `...`, with
  the same defaults as `dgp_birm_rs()` and `fit_stan()`. Misspelled or unknown
  arguments now raise an error; previously they were silently ignored and the
  default value was used instead.
* If no `seed` is given, `one_simulation()` now uses the same random seed for
  data generation and model fitting (previously two different random seeds).
  Results for a given `seed` are unchanged.

## Bug fixes

* `dgp_birm_rs()` works for a single latent trait (`theta_n = 1`), which
  previously failed.
* `sim_fun()` passes `adapt_delta` from `sim_grid` on to the model fit; in
  0.1.0 it was silently dropped and the default of 0.9 was used.
* `sim_fun()` no longer carries warnings over from one replication to later
  ones handled by the same worker.
* `sim_fun()` keeps the column types of `sim_grid` in the returned data frame.
* `sim_fun()` computes the number of workers for `workers = "half"` and
  `"quarter"` correctly.
* The example script (`inst/examples/sim_example.R`) assigned the `condition`
  numbers incorrectly; it now uses `rep(..., times = )` to match the row order
  of `expand.grid()`.

## Documentation

* Corrected outdated descriptions, among others for `summarize_fit()` (which
  reports correlation, bias, absolute bias and RMSE, not posterior SDs,
  minima and maxima) and for the initial values used by `fit_stan()`.
* Documented the required columns of `sim_grid` in `sim_fun()`.
* Fixed the examples of `one_simulation()` and `fit_stan()`.

## Package infrastructure

* The package now declares `Depends: R (>= 4.1.0)`, as it uses the native pipe
  `|>`.
* Functions from `stats` and `utils` are now imported explicitly;
  `R CMD check` passes without notes.
