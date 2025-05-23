#' Simulates data, fits a model and summarizes the fit
#'
#' Summarizes a fit object by computing means, standard deviations, minima and maxima of the posterior distributions of the model parameters, and the root mean squared error of the model parameter estimates.
#' @param ... arguments passed to \code{\link{dgp_birm_rs}}, \code{\link{fit_stan}} and \code{\link{summarize_fit}}
#' @export

one_simulation <- function(shorten = TRUE, ...){

  args <- list(...)

  dgp_birm_rs_args <- args[which(names(args) %in% c("n", "item_n", "theta_n", "var_thetas", "var_ers", "var_ars", "x_num", "cor_thetas", "cor_ers", "seed"))]
  fit_stan_args <- args[which(names(args) %in% c("stan_model", "theta_n", "seed", "iter", "chains", "adapt_delta", "warmup", "ars_prior", "init_vals"))]

  # write a list with the simulation conditions
  # conditions <- args

  data_list <- do.call(dgp_birm_rs, dgp_birm_rs_args)

  data <- data_list$df
  items <- data_list$items

  fit_stan_args$data <- data
  fit_stan_args$x_vec <- items$x
  fit_stan_args$T_vec <- args$item_n
  fit_stan_args$prefix <- "observed"


  # add tryCatch here?
  fit <- do.call(fit_stan, fit_stan_args)

  results <- summarize_fit(fit, data_list)
  rm(fit)

  # create a data.frame with one row and the args as columns. For args elements which are vectors, split them up
  # conditions <- data.frame(t(unlist(args)))
  # results$summary <- cbind(conditions, results$summary)


  results$data <- data
  results$items <- items

  # model_summary[["warnings"]] <- c(names(warnings_tmp), names(warnings()))
  # and here? or parent function.
  

  
  return(results)

}
