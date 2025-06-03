#' Simulates data, fits a model, and summarizes the fit
#'
#' This function acts as a wrapper for \code{\link{dgp_birm_rs}}, \code{\link{fit_stan}}, and \code{\link{summarize_fit}}. It first simulates data based on the specified parameters, fits the BIRM-RS model to the data, and finally summarizes the fit, including posterior summaries and RMSEs.
#'
#' @param n Number of observations (from \code{\link{dgp_birm_rs}}).
#' @param item_n Number of items or vector of items per trait; "auto" defaults to 10 (from \code{\link{dgp_birm_rs}}).
#' @param theta_n Number of latent traits (from \code{\link{dgp_birm_rs}} and \code{\link{fit_stan}}).
#' @param var_thetas Variances of latent traits (from \code{\link{dgp_birm_rs}}).
#' @param var_ers Variance of ERS parameter (from \code{\link{dgp_birm_rs}}).
#' @param var_ars Variance of ARS parameter (from \code{\link{dgp_birm_rs}}).
#' @param x_num Number of reverse-scored items per trait (from \code{\link{dgp_birm_rs}}).
#' @param cor_thetas Correlations between latent traits (from \code{\link{dgp_birm_rs}}).
#' @param cor_ers Correlations between ERS and latent traits (from \code{\link{dgp_birm_rs}}).
#' @param include_ERS Logical; whether to include ERS (from \code{\link{dgp_birm_rs}}).
#' @param include_ARS Logical; whether to include ARS (from \code{\link{dgp_birm_rs}}).
#' @param seed Random seed (used by all functions).
#'
#' @param stan_model Path to the Stan model file (from \code{\link{fit_stan}}).
#' @param x_vec Vector indicating reverse-scored items (from \code{\link{fit_stan}}).
#' @param T_vec Vector of item numbers per trait; "auto" divides total items by theta_n (from \code{\link{fit_stan}}).
#' @param iter Number of iterations excluding warmup (from \code{\link{fit_stan}}).
#' @param warmup Number of warmup iterations (from \code{\link{fit_stan}}).
#' @param chains Number of MCMC chains (from \code{\link{fit_stan}}).
#' @param adapt_delta Target acceptance rate for adaptation (from \code{\link{fit_stan}}).
#' @param prefix Prefix of variables to analyze; NULL uses all columns (from \code{\link{fit_stan}}).
#' @param ars_prior Prior for ARS parameter (from \code{\link{fit_stan}}).
#' @param init_vals Logical; whether to use initial values (from \code{\link{fit_stan}}).
#'
#' @param df_list List containing the data and items used to fit the model (from \code{\link{summarize_fit}}).
#' 
#' @param ... Additional arguments passed to \code{\link{dgp_birm_rs}}, \code{\link{fit_stan}}, and \code{\link{summarize_fit}}.
#' @export


one_simulation <- function(...){

  args <- list(...)

  # sort the arguments between the functions
  dgp_birm_rs_args <- args[which(names(args) %in% c("n", "item_n", "theta_n", "var_thetas", "var_ers", "var_ars", "x_num", "cor_thetas", "cor_ers", "seed"))]
  fit_stan_args <- args[which(names(args) %in% c("stan_model", "theta_n", "seed", "iter", "chains", "adapt_delta", "warmup", "ars_prior", "init_vals"))]

  # create the data
  data_list <- do.call(dgp_birm_rs, dgp_birm_rs_args)

  data <- data_list$df
  items <- data_list$items

  # add a few arguments based on the data
  fit_stan_args$data <- data
  fit_stan_args$x_vec <- items$x
  fit_stan_args$T_vec <- args$item_n
  fit_stan_args$prefix <- "observed"


  # fit the model
  fit <- do.call(fit_stan, fit_stan_args)

  # summarize the model 
  results <- summarize_fit(fit, data_list)
  rm(fit)

  # add the data and items to the results object (list)
  results$data <- data
  results$items <- items

  return(results)

}
