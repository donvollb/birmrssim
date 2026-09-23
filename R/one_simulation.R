#' Simulates data, fits a model, and summarizes the fit
#'
#' This function acts as a wrapper for \code{\link{dgp_birm_rs}}, \code{\link{fit_stan}}, and \code{\link{summarize_fit}}. It first simulates data based on the specified parameters, fits the BIRM-RS model to the data, and finally summarizes the fit, including posterior summaries and RMSEs.
#'
#' @inheritParams dgp_birm_rs
#' @inheritParams fit_stan
#' @param theta_n Number of latent traits (passed to \code{\link{dgp_birm_rs}} and \code{\link{fit_stan}}).
#' @param seed Seed for the random number generator (passed to \code{\link{dgp_birm_rs}} and \code{\link{fit_stan}}).
#'
#' @return A named list containing:
#' \describe{
#'   \item{summary}{Summary of simulation conditions, model fit and performance measures}
#'   \item{post_means}{Posterior means of parameters}
#'   \item{diagnostics}{MCMC diagnostics}
#'   \item{data}{The generated person parameters and simulated responses}
#'   \item{items}{The item parameters used in the simulation}
#' }
#' @details
#' The arguments \code{n} to \code{cor_ers} are passed to
#' \code{\link{dgp_birm_rs}}, the arguments \code{stan_model} to
#' \code{init_vals} to \code{\link{fit_stan}}, and \code{theta_n} and
#' \code{seed} to both. The defaults are the same as in these functions.
#'
#' The arguments \code{x_vec}, \code{T_vec}, and \code{prefix} (passed to
#' \code{\link{fit_stan}}) and \code{df_list} (passed to
#' \code{\link{summarize_fit}}) are set automatically from the simulated data
#' and cannot be supplied by the user.
#'
#'
#' @examples
#' \dontrun{
#' result <- one_simulation(
#'   n = 200, theta_n = 2, item_n = c(10, 10),
#'   stan_model = system.file("stan", "BIRM_RS.stan", package = "birmrssim"),
#'   iter = 1000, warmup = 500, chains = 2
#' )
#' }
#'
#' @seealso \code{\link{dgp_birm_rs}}, \code{\link{fit_stan}}, \code{\link{summarize_fit}}
#'
#' @export


one_simulation <- function(n = 2000, item_n = "auto", theta_n = 1,
  var_thetas = "auto", var_ers = 0.3, var_ars = 0.3, cor_thetas = "auto",
  x_num = "auto", cor_ers = "auto", seed = sample(1:1e9, 1),
  stan_model, iter = 8000, warmup = 4000, chains = 6, adapt_delta = 0.9,
  ars_prior = 0.9, init_vals = TRUE) {

  # create the data
  data_list <- dgp_birm_rs(
    n = n, item_n = item_n, theta_n = theta_n, var_thetas = var_thetas,
    var_ers = var_ers, var_ars = var_ars, cor_thetas = cor_thetas,
    x_num = x_num, cor_ers = cor_ers, seed = seed
  )

  data <- data_list$df
  items <- data_list$items

  # fit the model; x_vec, T_vec and prefix follow from the simulated data
  # (item_n = "auto" gives equally many items per trait, so T_vec = "auto" matches)
  fit <- fit_stan(
    data, stan_model = stan_model, theta_n = theta_n, x_vec = items$x,
    T_vec = item_n, iter = iter, warmup = warmup, chains = chains, seed = seed,
    adapt_delta = adapt_delta, prefix = "observed", ars_prior = ars_prior,
    init_vals = init_vals
  )

  # summarize the model
  results <- summarize_fit(fit, data_list)
  rm(fit)

  # add the data and items to the results object (list)
  results$data <- data
  results$items <- items

  return(results)

}
