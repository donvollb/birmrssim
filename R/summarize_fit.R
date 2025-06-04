#' Summarize a fit object from \code{\link{fit_stan}}
#'
#' Summarizes the results from a fit object, including means, standard deviations, minima, and maxima of the posterior distributions of the model parameters, as well as the root mean squared error (RMSE) of the model parameter estimates.
#' @param fit Fit object from \code{\link{fit_stan}}.
#' @param df_list List containing the data used to fit the model and the items used to fit the model.
#' @param theta_n Number of latent traits. If not specified, it will be automatically inferred from the column names of \code{df_list$df}.
#' @return A list containing summary statistics and the RMSEs of the model parameters.
#' @import loo
#' @import dplyr
#' @export


summarize_fit <- function(fit, df_list, theta_n = "auto") {

  # get true values
  true <- df_list$df
  items <- df_list$items

  # get the number of thetas if not specified
  if (theta_n == "auto") grepl("theta", colnames(true)) |> sum() -> theta_n

  # Custom RMSE function
  rmse <- function(post_val, true_val) {
    sqrt(mean((post_val - true_val)^2))
  }

  # Extract relevant posterior means
  means <- fit$summary(variables = c("theta", "ars", "delta", "tau", "L_Sigma")) |>
    select(variable, mean)

  # Dynamically extract thetas
  theta_estimates <- list()
  for (i in 1:theta_n) {
    var_name <- paste0("theta_", i)
    theta_estimates[[var_name]] <- means[grepl(paste0("theta.*," , i, "]"), means$variable), "mean"] |> unlist()
  }

  # ERS is the last column of theta
  ers <- means[grepl(paste0("theta.*," , theta_n + 1, "]"), means$variable), "mean"] |> unlist()
  ars <- means[grepl("ars", means$variable), "mean"] |> unlist()
  delta <- means[grepl("delta", means$variable), "mean"] |> unlist()
  tau <- means[grepl("tau", means$variable), "mean"] |> unlist()

  # Start results list
  result <- data.frame(matrix(nrow = 1, ncol = 0))

  # Add correlations, biases, variances, means, rmse for theta
  for (i in 1:theta_n) {
    name <- paste0("theta", i)
    est <- theta_estimates[[paste0("theta_", i)]]
    tru <- true[[name]]

    result[[paste0("est_cor_", name)]] <- cor(est, tru)
    result[[paste0("est_bias_", name)]] <- mean(est - tru)
    result[[paste0("est_bias_abs_", name)]] <- abs(est - tru) |> mean()
    result[[paste0("est_var_", name)]] <- var(est)
    result[[paste0("est_mean_", name)]] <- mean(est)
    result[[paste0("est_md_", name)]] <- median(est)
    result[[paste0("est_rmse_", name)]] <- rmse(est, tru)
  }

  # Add stats for ERS
  result$est_cor_ers <- cor(ers, true$ers)
  result$est_bias_ers <- mean(ers - true$ers)
  result$est_bias_abs_ers <- abs(ers - true$ers) |> mean()
  result$est_var_ers <- var(ers)
  result$est_mean_ers <- mean(ers)
  result$est_md_ers <- median(ers)
  result$est_rmse_ers <- rmse(ers, true$ers)

  # Add stats for ARS
  result$est_cor_ars <- cor(ars, true$ars)
  result$est_bias_ars <- mean(ars - true$ars)
  result$est_bias_abs_ars <- abs(ars - true$ars) |> mean()
  result$est_var_ars <- var(ars)
  result$est_mean_ars <- mean(ars)
  result$est_md_ars <- median(ars)
  result$est_rmse_ars <- rmse(ars, true$ars)

  # Add correlations between person parameters (important: these are the correlations based on the posterior means)
  param_df <- data.frame(
    ers = ers,
    ars = ars
  )
  for (i in 1:theta_n) {
    param_df[[paste0("theta", i)]] <- theta_estimates[[paste0("theta_", i)]]
  }

  param_names <- names(param_df)
  for (i in 1:(length(param_names) - 1)) {
    for (j in (i + 1):length(param_names)) {
      name_i <- param_names[i]
      name_j <- param_names[j]
      result[[paste0("est_cor_", name_i, "_", name_j)]] <- cor(param_df[[name_i]], param_df[[name_j]])
    }
  }

  # now add the estimated variance-covariance matrix for theta_n and ERS
  # the covariances in one diagonal block are 0 by default
  # we do not remove them here to keep compability with different numbers of latent traits
  covs <- means[grepl("L_Sigma", means$variable), ]
  covs_row <- t(covs$mean)
  colnames(covs_row) <- covs$variable

  # add it to result
  result <- cbind(result, covs_row)

  # Add stats for item parameters
  result$est_cor_delta <- cor(delta, items$delta)
  result$est_bias_delta <- mean(delta - items$delta)
  result$est_bias_abs_delta <- abs(delta - items$delta) |> mean()
  result$est_rmse_delta <- rmse(delta, items$delta)

  result$est_cor_tau <- cor(tau, items$tau)
  result$est_bias_tau <- mean(tau - items$tau)
  result$est_bias_abs_tau <- abs(tau - items$tau) |> mean()
  result$est_rmse_tau <- rmse(tau, items$tau)

  # Compute LOO and WAIC
  log_lik_matrix <- fit$draws("log_lik", format = "draws_matrix")
  loo_result <- loo(log_lik_matrix)
  waic_result <- waic(log_lik_matrix)

  result$elpd_loo <- loo_result$estimates["elpd_loo", "Estimate"]
  result$p_loo <- loo_result$estimates["p_loo", "Estimate"]
  result$looic <- loo_result$estimates["looic", "Estimate"]
  result$elpd_waic <- waic_result$estimates["elpd_waic", "Estimate"]
  result$p_waic <- waic_result$estimates["p_waic", "Estimate"]
  result$waic <- waic_result$estimates["waic", "Estimate"]

  # Compute proportion of problematic Pareto k values (k > 0.7)
  k_values <- loo_result$diagnostics$pareto_k
  result$pareto_k_bad_prop <- mean(k_values > 0.7, na.rm = TRUE)

  
  waic_p <- waic_result$pointwise[, "p_waic"]
  result$p_waic_bad_prop <- mean(waic_p > 0.4, na.rm = TRUE)

  # Compute ESS and Rhat
  ess_hats <- fit$summary(variables = NULL, c("rhat", "ess_bulk", "ess_tail"))

  # Return both result row and posterior means
  return(list(
    summary = result,
    post_means = c(theta_estimates, list(ers = ers, ars = ars, delta = delta, tau = tau)),
    diagnostics = ess_hats
  ))
}
