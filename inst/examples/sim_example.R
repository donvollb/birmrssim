# =============================================================================
# Example Simulation Script for the birmrssim Package
# =============================================================================
# This script demonstrates how to use the birmrssim package to run a small
# simulation study evaluating parameter recovery in the BIRM-RS model.
# For a full simulation study, increase the number of conditions and replications.
# =============================================================================

# Load packages
library(cmdstanr)
library(dplyr)
library(loo)
library(future)
library(future.apply)
library(progressr)
library(birmrssim)

# -----------------------------------------------------------------------------
# Step 1: Access the Stan model files included in the package
# -----------------------------------------------------------------------------
model_string <- system.file("stan", "BIRM_RS.stan", package = "birmrssim")

# -----------------------------------------------------------------------------
# Step 2: Define simulation conditions
# -----------------------------------------------------------------------------
conditions_list <- list(
  n = c(300, 600),            # Sample size
  theta_n = 2,                # Number of latent traits
  item_n1 = 6,                # Number of items for first latent trait
  item_n2 = 6,                # Number of items for second latent trait
  x_num1 = c(0, 3),           # Number of reverse-scored items, first trait
  x_num2 = c(0, 3),           # Number of reverse-scored items, second trait
  var_thetas1 = 1,            # Variance of first latent trait
  var_thetas2 = 1,            # Variance of second latent trait
  var_ers = c(0.5, 0.3),      # Variance of ERS parameter
  var_ars = c(0.5, 0.3),      # Variance of ARS parameter
  cor_ers1 = 0.1,             # Correlation of ERS and first latent trait
  cor_ers2 = 0.1,             # Correlation of ERS and second latent trait
  cor_thetas = 0.2,           # Correlation between latent traits
  ars_prior = 0.9,            # Prior SD for ARS parameter
  chains = 2,                 # Number of MCMC chains
  iter = 1000,                # Post-warmup iterations (increase for real study)
  warmup = 500,               # Warmup iterations
  adapt_delta = 0.9,          # Adaptation delta
  init_vals = FALSE,          # Use initial values?
  stan_model = model_string,  # Path to Stan model
  seed = NA,                  # Seeds will be added below
  replication = 1:3           # Number of replications (increase for real study)
)

# -----------------------------------------------------------------------------
# Step 3: Build the simulation grid
# -----------------------------------------------------------------------------
sim_grid <- expand.grid(conditions_list)
sim_grid$stan_model <- as.character(sim_grid$stan_model)

# Remove conditions where var_ers != var_ars
sim_grid <- sim_grid[sim_grid$var_ers == sim_grid$var_ars, ]

# Add condition and index columns
sim_grid$condition <- rep(1:(nrow(sim_grid) / length(conditions_list$replication)),
                          each = length(conditions_list$replication))
sim_grid$index <- 1:nrow(sim_grid)
sim_grid <- sim_grid[, c("index", "condition", 
                          setdiff(names(sim_grid), c("index", "condition")))]

# Add unique seeds
sim_grid$seed <- sample(1:1e7, size = nrow(sim_grid), replace = FALSE)

# -----------------------------------------------------------------------------
# Step 4: Run the simulations
# -----------------------------------------------------------------------------
# Use parallel_type = "multisession" on Windows or Linux,
# "multicore" on Linux/macOS only, or "sequential" for testing.
results_df <- sim_fun(
  sim_grid      = sim_grid,
  results_path  = "results",   # Directory to save individual result files, has to exist before running
  parallel_type = "sequential" # Change to "multisession" for parallel execution
)

# -----------------------------------------------------------------------------
# Step 5: Save and inspect results
# -----------------------------------------------------------------------------
save(results_df, file = "results/results_example.rda")
print(head(results_df))
print(sessionInfo())