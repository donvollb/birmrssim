#' Fit a BIRM-RS model
#'
#' Fits a BIRM-RS model to the provided dataset.
#'
#' @param data Data frame containing the observed data.
#' @param stan_model Path to the Stan model file.
#' @param theta_n Number of latent traits.
#' @param x_vec Vector indicating reverse-scored items (1 = normally scored, -1 = reverse-scored).
#' @param T_vec Vector specifying the number of items per trait. Defaults to "auto", which sets it to the total number of items divided by theta_n.
#' @param iter Number of iterations (excluding warmup).
#' @param warmup Number of warmup iterations.
#' @param chains Number of MCMC chains.
#' @param seed Random seed for reproducibility.
#' @param adapt_delta Target acceptance probability for step size adaptation in Stan. Defaults to 0.9.
#' @param prefix Optional prefix indicating which variables to analyze. If NULL, all columns are used. Default is "observed".
#' @param ars_prior Prior for the ARS parameter. Default is 0.9.
#' @param init_vals Logical; indicates whether to use initial values. Default is TRUE.
#'
#' @details
#' Response values are clamped to the interval \[0.0001, 0.9999\] before fitting
#' to avoid boundary issues, following the approach of Noel & Dauvier (2007).
#' 
#' The Stan model files included in this package can be accessed via
#' \code{system.file()}. Available models are:
#' \describe{
#'   \item{BIRM_RS.stan}{Full model with ERS and ARS}
#'   \item{ERS_only.stan}{Model with ERS only}
#'   \item{ARS_only.stan}{Model with ARS only}
#'   \item{no_RS.stan}{Model without response styles}
#' }
#' Example: \code{system.file("stan", "BIRM_RS.stan", package = "birmrssim")}
#'
#' If \code{init_vals = TRUE}, initial values are drawn randomly for all chain
#' parameters to facilitate convergence.
#'
#' @return A \code{CmdStanMCMC} object (from \pkg{cmdstanr}) representing the fitted model.
#' 
#' @references
#' Noel, Y., & Dauvier, B. (2007). A Beta Item Response Model for Continuous
#' Bounded Responses. \emph{Applied Psychological Measurement}, \emph{31}(1),
#' 47--73. \doi{10.1177/0146621605287691}
#'
#' @import cmdstanr
#' @examples
#' \dontrun{
#' model_path <- system.file("stan", "BIRM_RS.stan", package = "birmrssim")
#' results <- fit_stan(data, stan_model = model_path, theta_n = 2,
#'                     x_vec = c(rep(1, 4), rep(-1, 4)), iter = 8000, warmup = 4000, chains = 6,
#'                     seed = sample(1:1e9, 1), adapt_delta = 0.9, prefix = "observed", init_vals = TRUE)
#' }
#' @export


fit_stan <- function(
    data, stan_model, theta_n, x_vec, T_vec = "auto", 
    iter = 8000, warmup = 4000, chains = 6, seed = sample(1:1e9, 1), adapt_delta = 0.9, 
    prefix = "observed", ars_prior = 0.9, init_vals = TRUE) {

  
  # if there is no prefix, use all columns
  if (!is.null(prefix)) {data_cols <- grepl(prefix, colnames(data))} else {data_cols <- rep(TRUE, ncol(data))}
  
  data <- data[, data_cols]
  
  # if T_vec is auto, set it to the number of columns divided by theta_n
  # that would mean that each theta would have the same number of items
  if (T_vec[1] == "auto") T_vec <- rep(ncol(data)/theta_n, theta_n)
  
  

  data[data > .9999] <- 0.9999
  data[data < .0001] <- 0.0001
  # see paper from Noel & Dauvier

  n <- nrow(data)

  # Package into a list for Stan
  if (theta_n == 1) stan_data <- list(N = n, T = array(T_vec, dim = c(1)), x = x_vec, r = as.matrix(data), theta_n = theta_n, ars_prior = ars_prior)
  if (theta_n > 1) stan_data <- list(N = n, T = T_vec, x = x_vec, r = as.matrix(data), theta_n = theta_n, ars_prior = ars_prior)

  set.seed(seed)

  # load the model (and compile if necessary)
  stan_model_path <- stan_model
  stan_model <- cmdstan_model(stan_model_path)

  # initial values (if relevant)
  if (init_vals == TRUE) {

    # detect ARS from compiled model parameter names
    has_ars <- "ars" %in% names(stan_model$variables()$parameters)
    # detect ERS from Stan source (ERS models use theta_n + 1 as theta dimension)
    has_ers    <- any(grepl("theta_n\\s*\\+\\s*1", readLines(stan_model_path, warn = FALSE)))
    theta_cols <- if (has_ers) theta_n + 1 else theta_n

    common_inits <- list(
      theta      = matrix(rnorm(n * theta_cols, 0, 1), nrow = n, ncol = theta_cols),
      delta      = rnorm(ncol(data), 0, 5),
      tau        = rep(1.5, ncol(data)),
      sigma      = rep(1, theta_cols),
      Sigma_corr = diag(theta_cols)
    )
    if (has_ars) common_inits$ars <- rnorm(n, 0, 1)

    init <- rep(list(common_inits), chains)

  } else {
    init <- NULL
  }


  # fit the stan model
  fit <- stan_model$sample(
    data = stan_data,
    seed = seed,
    chains = chains,
    iter_warmup = warmup,
    iter_sampling = iter,
    adapt_delta = adapt_delta,
    parallel_chains = chains,
    init = init,
  )

  return(fit)
}

