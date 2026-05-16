#' Summarize a fit object from \code{\link{fit_stan}}
#'
#' Summarizes the results from a fit object, including means, standard deviations,
#' minima, and maxima of the posterior distributions of the model parameters, as
#' well as the root mean squared error (RMSE) of the model parameter estimates.
#'
#' @param fit Fit object from \code{\link{fit_stan}}.
#' @param df_list List containing the data and item parameters used to fit the model,
#'   as returned by \code{\link{dgp_birm_rs}}.
#' @param theta_n Number of latent traits. If not specified, it will be automatically
#'   inferred from the column names of \code{df_list$df}.
#' @details
#' The presence of ERS and ARS is detected automatically from the variable names
#' in the fitted model: ERS is assumed present when \code{theta} has
#' \code{theta_n + 1} columns (last column = ERS), and ARS is assumed present
#' when an \code{ars} parameter exists. Columns for absent RS parameters are
#' set to \code{NA}, and the \code{L_Sigma} block is always output as a
#' \code{(theta_n+1) x (theta_n+1)} matrix (ERS row/column padded with \code{NA}
#' when ERS is absent). This keeps the output structure identical across all four
#' model variants so that results can be combined with \code{dplyr::bind_rows()}.
#'
#' Performance measures (correlations, bias, RMSE) are computed by comparing
#' posterior means to the true values used in the data generating process.
#' Correlations between person parameters (ERS, ARS, latent traits) are based
#' on posterior means rather than full posterior distributions.
#'
#' Model fit is assessed via LOO-CV and WAIC using the \pkg{loo} package.
#' The proportion of problematic Pareto-k values (k > 0.7) and problematic
#' p_waic values (> 0.4) are also reported.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{summary}{A one-row data frame containing performance measures
#'   (correlations, bias, RMSE) for all person and item parameters, estimated
#'   correlations among person parameters, Cholesky factor elements of the
#'   covariance matrix, and model fit indices (LOO-CV, WAIC).}
#'   \item{post_means}{A named list of posterior means for all person and item
#'   parameters (latent traits, delta, tau, and ERS/ARS when present in the
#'   fitted model).}
#'   \item{diagnostics}{A data frame containing Rhat, bulk ESS, and tail ESS
#'   for all model parameters.}
#' }
#'
#' @import loo
#' @import dplyr
#' @export


summarize_fit <- function(fit, df_list, theta_n = "auto") {

  # get true values
  true <- df_list$df
  items <- df_list$items

  # get the number of thetas if not specified
  if (theta_n == "auto") grepl("theta", colnames(true)) |> sum() -> theta_n

  # Auto-detect whether ERS and ARS are present in the fitted model
  all_vars <- dimnames(fit$draws())[[3]]
  include_ARS <- any(grepl("^ars\\[", all_vars))
  include_ERS <- any(grepl(paste0("^theta\\[.*,", theta_n + 1, "\\]"), all_vars))

  # Custom RMSE function
  rmse <- function(post_val, true_val) {
    sqrt(mean((post_val - true_val)^2))
  }

  # Build variable list based on which RS are included
  vars <- c("theta", "delta", "tau", "L_Sigma")
  if (include_ARS) vars <- c(vars, "ars")

  # Extract relevant posterior means
  means <- fit$summary(variables = vars) |>
    select(variable, mean)

  # Dynamically extract thetas
  theta_estimates <- list()
  for (i in 1:theta_n) {
    var_name <- paste0("theta_", i)
    theta_estimates[[var_name]] <- means[grepl(paste0("theta.*,", i, "]"), means$variable), "mean"] |> unlist()
  }

  # ERS is the last column of theta (theta_n + 1)
  if (include_ERS) {
    ers <- means[grepl(paste0("theta.*,", theta_n + 1, "]"), means$variable), "mean"] |> unlist()
  }
  if (include_ARS) {
    ars <- means[grepl("ars", means$variable), "mean"] |> unlist()
  }

  delta <- means[grepl("delta", means$variable), "mean"] |> unlist()
  tau   <- means[grepl("tau",   means$variable), "mean"] |> unlist()

  # Start results list
  result <- data.frame(matrix(nrow = 1, ncol = 0))

  # Add correlations, biases, variances, means, rmse for theta
  for (i in 1:theta_n) {
    name <- paste0("theta", i)
    est <- theta_estimates[[paste0("theta_", i)]]
    tru <- true[[name]]

    result[[paste0("est_cor_", name)]]      <- cor(est, tru)
    result[[paste0("est_bias_", name)]]     <- mean(est - tru)
    result[[paste0("est_bias_abs_", name)]] <- abs(est - tru) |> mean()
    result[[paste0("est_var_", name)]]      <- var(est)
    result[[paste0("est_mean_", name)]]     <- mean(est)
    result[[paste0("est_md_", name)]]       <- median(est)
    result[[paste0("est_rmse_", name)]]     <- rmse(est, tru)
  }

  # Add stats for ERS
  if (include_ERS) {
    result$est_cor_ers      <- cor(ers, true$ers)
    result$est_bias_ers     <- mean(ers - true$ers)
    result$est_bias_abs_ers <- abs(ers - true$ers) |> mean()
    result$est_var_ers      <- var(ers)
    result$est_mean_ers     <- mean(ers)
    result$est_md_ers       <- median(ers)
    result$est_rmse_ers     <- rmse(ers, true$ers)
  } else {
    result$est_cor_ers      <- NA
    result$est_bias_ers     <- NA
    result$est_bias_abs_ers <- NA
    result$est_var_ers      <- NA
    result$est_mean_ers     <- NA
    result$est_md_ers       <- NA
    result$est_rmse_ers     <- NA
  }

  # Add stats for ARS
  if (include_ARS) {
    result$est_cor_ars      <- cor(ars, true$ars)
    result$est_bias_ars     <- mean(ars - true$ars)
    result$est_bias_abs_ars <- abs(ars - true$ars) |> mean()
    result$est_var_ars      <- var(ars)
    result$est_mean_ars     <- mean(ars)
    result$est_md_ars       <- median(ars)
    result$est_rmse_ars     <- rmse(ars, true$ars)
  } else {
    result$est_cor_ars      <- NA
    result$est_bias_ars     <- NA
    result$est_bias_abs_ars <- NA
    result$est_var_ars      <- NA
    result$est_mean_ars     <- NA
    result$est_md_ars       <- NA
    result$est_rmse_ars     <- NA
  }

  # Add correlations between person parameters
  # Build a helper data frame; absent parameters become NA columns so the
  # correlation loop always produces a column (with NA as value)
  param_df <- data.frame(
    ers = if (include_ERS) ers else rep(NA_real_, nrow(true)),
    ars = if (include_ARS) ars else rep(NA_real_, nrow(true))
  )
  for (i in 1:theta_n) {
    param_df[[paste0("theta", i)]] <- theta_estimates[[paste0("theta_", i)]]
  }

  param_names <- names(param_df)
  for (i in 1:(length(param_names) - 1)) {
    for (j in (i + 1):length(param_names)) {
      name_i <- param_names[i]
      name_j <- param_names[j]
      val <- if (all(!is.na(param_df[[name_i]])) && all(!is.na(param_df[[name_j]]))) {
        cor(param_df[[name_i]], param_df[[name_j]])
      } else {
        NA_real_
      }
      result[[paste0("est_cor_", name_i, "_", name_j)]] <- val
    }
  }

  # L_Sigma: always output a (theta_n+1) x (theta_n+1) block.
  # When ERS is absent the model only estimates a theta_n x theta_n matrix;
  # the entries for the ERS dimension (row/col theta_n+1) are padded with NA.
  covs <- means[grepl("L_Sigma", means$variable), ]

  if (include_ERS) {
    covs_row <- t(covs$mean)
    colnames(covs_row) <- covs$variable
  } else {
    full_dim   <- theta_n + 1
    full_names <- paste0("L_Sigma[", rep(1:full_dim, times = full_dim), ",",
                         rep(1:full_dim, each  = full_dim), "]")
    covs_full  <- setNames(rep(NA_real_, length(full_names)), full_names)
    covs_full[covs$variable] <- covs$mean
    covs_row   <- t(covs_full)
  }

  result <- cbind(result, covs_row)

  # Add stats for item parameters
  result$est_cor_delta      <- cor(delta, items$delta)
  result$est_bias_delta     <- mean(delta - items$delta)
  result$est_bias_abs_delta <- abs(delta - items$delta) |> mean()
  result$est_rmse_delta     <- rmse(delta, items$delta)

  result$est_cor_tau      <- cor(tau, items$tau)
  result$est_bias_tau     <- mean(tau - items$tau)
  result$est_bias_abs_tau <- abs(tau - items$tau) |> mean()
  result$est_rmse_tau     <- rmse(tau, items$tau)

  # Compute LOO and WAIC
  log_lik_matrix <- fit$draws("log_lik", format = "draws_matrix")
  loo_result  <- loo(log_lik_matrix)
  waic_result <- waic(log_lik_matrix)

  result$elpd_loo  <- loo_result$estimates["elpd_loo",  "Estimate"]
  result$p_loo     <- loo_result$estimates["p_loo",     "Estimate"]
  result$looic     <- loo_result$estimates["looic",     "Estimate"]
  result$elpd_waic <- waic_result$estimates["elpd_waic", "Estimate"]
  result$p_waic    <- waic_result$estimates["p_waic",    "Estimate"]
  result$waic      <- waic_result$estimates["waic",      "Estimate"]

  k_values <- loo_result$diagnostics$pareto_k
  result$pareto_k_bad_prop <- mean(k_values > 0.7, na.rm = TRUE)

  waic_p <- waic_result$pointwise[, "p_waic"]
  result$p_waic_bad_prop <- mean(waic_p > 0.4, na.rm = TRUE)

  # Compute ESS and Rhat
  ess_hats <- fit$summary(variables = NULL, c("rhat", "ess_bulk", "ess_tail"))

  # Build post_means; omit absent RS parameters
  post_means <- theta_estimates
  if (include_ERS) post_means$ers <- ers
  if (include_ARS) post_means$ars <- ars
  post_means$delta <- delta
  post_means$tau   <- tau

  return(list(
    summary    = result,
    post_means = post_means,
    diagnostics = ess_hats
  ))
}
