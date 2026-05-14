# birmrssim

An R package for evaluating parameter recovery of the **Beta Item Response Model with Response Styles (BIRM-RS)**. The package provides functions to simulate data from the BIRM-RS, fit the model via Stan, and summarize parameter recovery across customizable simulation conditions. Stan model files for four model variants are included.

This package was developed as a supplement to the following paper:

> *Citation will be added upon publication.*

---

## Installation

The package can be installed directly from GitHub using `pak`:

```r
# install.packages("pak")
pak::pak("donvollb/birmrssim")
```

### Dependencies

The package relies on [`cmdstanr`](https://mc-stan.org/cmdstanr/) to fit Stan models. `cmdstanr` is not on CRAN and must be installed separately:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
```

After installing `cmdstanr`, install CmdStan itself:

```r
cmdstanr::install_cmdstan()
```

---

## Included Stan Models

Four Stan model files are included in the package and can be accessed via `system.file()`:

| File | Description |
|---|---|
| `BIRM_RS.stan` | Full model with ERS and ARS |
| `ERS_only.stan` | Model with ERS only |
| `ARS_only.stan` | Model with ARS only |
| `no_RS.stan` | Model without response styles |

```r
# Access a model file
model_path <- system.file("stan", "BIRM_RS.stan", package = "birmrssim")
```

---

## Overview of Main Functions

| Function | Description |
|---|---|
| `dgp_birm_rs()` | Simulate data from the BIRM-RS |
| `fit_stan()` | Fit a BIRM-RS model via Stan |
| `summarize_fit()` | Summarize parameter recovery from a fitted model |
| `one_simulation()` | Wrapper: simulate, fit, and summarize in one step |
| `sim_fun()` | Run simulations in parallel across a condition grid |

---

## Quick Start

```r
library(birmrssim)

# Access the Stan model included in the package
model_path <- system.file("stan", "BIRM_RS.stan", package = "birmrssim")

# Run a single simulation
result <- one_simulation(
  n          = 300,
  theta_n    = 2,
  item_n     = c(6, 6),
  stan_model = model_path,
  iter       = 1000,
  warmup     = 500,
  chains     = 2,
  seed       = 42
)

# Inspect results
result$summary
result$diagnostics
```

---

## Example Simulation Study

A fully worked example script is included in the package:

```r
# View the example script
example_path <- system.file("examples", "sim_example.R", package = "birmrssim")
file.show(example_path)
```

The script demonstrates how to define a simulation grid, run simulations in parallel, and save results.

---

## Citation

If you use this package, please cite the associated paper:

> *Citation will be added upon publication.*

---

## License

GPL (>= 3)