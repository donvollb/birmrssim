#' Generate a dataset with the BIRM-RS
#'
#' Generates item and person parameters, then creates a dataset based on the BIRM-RS.
#' @param n Number of observations.
#' @param item_n Number of items. For multiple traits, this can be a vector specifying the number of items per trait. When set to "auto" (default), it defaults to 10 items per trait.
#' @param theta_n Number of latent traits (tested up to 3).
#' @param var_thetas Variances of the latent traits. Defaults to 1 ("auto") for all traits.
#' @param var_ers Variance of the ERS parameter. Defaults to 1.
#' @param var_ars Variance of the ARS parameter. Defaults to 1.
#' @param x_num Number of reverse-scored items per trait. Defaults to floor(item_n/2) ("auto").
#' @param cor_thetas Correlations among latent traits. For three traits, the order is: cor(theta1, theta2), cor(theta1, theta3), cor(theta2, theta3). Defaults to runif((theta_n*(theta_n - 1))/2, -0.4, 0.4) ("auto").
#' @param cor_ers Correlations between ERS and latent traits. For three traits, the order is: cor(theta1, ERS), cor(theta2, ERS), cor(theta3, ERS). Defaults to runif(theta_n, -0.3, 0.3).
#' @param include_ERS Logical; whether to include ERS in the model. Defaults to TRUE.
#' @param include_ARS Logical; whether to include ARS in the model. Defaults to TRUE.
#' @param seed Optional; seed for the random number generator.
#' @details This function generates item and person parameters based on the BIRM-RS, and then simulates a dataset from these parameters.
#' @return A list with two data frames: one containing item parameters, and the other containing person parameters and the simulated data.
#' @importFrom MASS mvrnorm
#' @examples

#' dfs_list <- dgp_birm_rs(n = 2000, item_n = 10, theta_n = 1, var_thetas = rep(1, theta_n), var_ers = 1, var_ars = 1, x_num = floor(item_n/2), cor_ers = 0.2)
#' @export

dgp_birm_rs <- function(n = 2000, item_n = "auto", theta_n = 1, 
  var_thetas = "auto", var_ers = 1, var_ars = 1, cor_thetas = "auto",
  x_num = "auto", cor_ers = "auto", include_ERS = TRUE, include_ARS = TRUE, 
  seed = sample(1:1e9, 1)) {

  set.seed(seed)

  # simulate values if none were set
  if (item_n[1] == "auto") item_n <- rep(10, theta_n)
  if (var_thetas[1] == "auto") var_thetas <- rep(1, theta_n)
  if (cor_thetas[1] == "auto") cor_thetas <- runif((theta_n*(theta_n - 1))/2, -0.4, 0.4)
  if (x_num[1] == "auto") x_num <- floor(item_n/2)
  if (cor_ers[1] == "auto") cor_ers <- runif(theta_n, -0.3, 0.3)
  # df <- data.frame(id = 1:n)

  # check if length(item_n) is equal to theta_n
  if (length(item_n) != theta_n) {
    stop("length(item_n) must be equal to theta_n")
  }

  # if length(x_num) is not equal to length(item_n), throw an error
  if (length(x_num) != length(item_n)) stop("length(x_num) must be equal to length(item_n)")

  # if length(var_thetas) is not equal to theta_n, throw an error
  if (length(var_thetas) != theta_n) stop("length(var_thetas) must be equal to theta_n")

  # if length(cor_thetas) is not equal to (theta_n*(theta_n - 1))/2, throw an error
  if (length(cor_thetas) != (theta_n*(theta_n - 1))/2) stop("length(cor_thetas) must be equal to (theta_n*(theta_n - 1))/2")

  # if length(cor_ers) is not equal to theta_n, throw an error
  if (length(cor_ers) != theta_n) stop("length(cor_ers) must be equal to theta_n")


  # 1. Create the covariance matrix for the latent traits (theta)
  CovTheta <- diag(var_thetas)  # Start with a diagonal matrix using the variances
  counter <- 1
  for (i in 1:(theta_n - 1)) {
    for (j in (i + 1):theta_n) {
      # Fill the off-diagonals with the appropriate covariance computed from cor_thetas
      CovTheta[i, j] <- cor_thetas[counter] * sqrt(var_thetas[i] * var_thetas[j])
      CovTheta[j, i] <- CovTheta[i, j]  # ensure symmetry
      counter <- counter + 1
    }
  }
  
  # 2. Compute the covariance between each theta and ERS
  # This creates a vector of covariances for each latent trait and ERS.
  cov_ers_vec <- cor_ers * sqrt(var_ers) * sqrt(var_thetas)
  
  # 3. Build the full covariance matrix (Sigma) for theta and ERS.
  # The matrix will be of dimension (theta_n + 1) x (theta_n + 1)
  Sigma <- matrix(0, nrow = theta_n + 1, ncol = theta_n + 1)
  Sigma[1:theta_n, 1:theta_n] <- CovTheta  # assign the latent traits covariance block
  Sigma[theta_n + 1, theta_n + 1] <- var_ers  # assign the ERS variance
  Sigma[1:theta_n, theta_n + 1] <- cov_ers_vec  # fill covariances between theta and ERS
  Sigma[theta_n + 1, 1:theta_n] <- cov_ers_vec  # ensure symmetry

  # start the dataframe
  df <- data.frame(
    id = 1:n
  )

  # simulate latent traits and ERS from the covariance matrix
  df[, c(paste0("theta", 1:theta_n), "ers")] <- MASS::mvrnorm(
    n = n,
    mu = rep(0, theta_n + 1),
    Sigma = Sigma
  )

  # simulate ARS
  df$ars <- rnorm(n, 0, sqrt(var_ars))


  # create variable for items
  items_total <- sum(item_n)
  
  # create x vector
  x_vec <- NULL
  for (i in 1:theta_n) {

    # create a vector that starts with 1 to item_n[i] - x_num and then -1 to item_n[i]
    tmp_vec <- rep(1, item_n[i] - x_num[i])
    tmp_vec <- c(tmp_vec, rep(-1, x_num[i]))
    x_vec <- c(x_vec, tmp_vec)     

  }

  # create item parameters
  items <- data.frame(
    item_id = 1:items_total,
    delta = c(0, runif(items_total - 1, -3, 3)),
    tau = runif(items_total, 0, 3),
    x = x_vec
  )

  # simulate item responses
  # helper
  e <- 0
  
  for (i in 1:theta_n) {
    col_name <- paste0("theta", i)
    
    for (k in 1:item_n[i]) {

      # which item?
      f <- e + k

      
      if (include_ARS) {
        middle <- df[, col_name] - items$delta[f] + items$x[f] * df$ars
      } else {
        middle <- df[, col_name] - items$delta[f]
      }
      
      if (include_ERS) {
        numerator_m <- exp(df$ers)*middle + items$tau[f]
        numerator_n <- -exp(df$ers)*middle + items$tau[f]
      } else {
        numerator_m <- middle + items$tau[f]
        numerator_n <- -middle + items$tau[f]
      }

      mj <- exp(numerator_m/2)
      nj <- exp(numerator_n/2)

      df[paste0("observed_", i, "_", k)] <- rbeta(n, mj, nj) 

    }
    
    e <- e + item_n[i]
  }

  # return data and item parameters
  return(list(df = df, items = items))

} 
