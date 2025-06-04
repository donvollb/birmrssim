#' Run Simulations in Parallel with Progress Reporting
#'
#' This function runs simulations in parallel based on the specified simulation grid.
#' It compiles Stan models, executes simulations across multiple workers, handles warnings and errors gracefully, and saves individual simulation results.
#'
#' @param sim_grid A data frame containing the grid of simulation parameters.
#' @param results_path File path to save individual simulation result files.
#' @param workers Number of parallel workers to use. Can be a numeric value or one of "half" or "quarter" to specify relative to available cores.
#' @param save_all Logical; whether to save all individual simulation result objects to disk.
#' @param parallel_type Type of parallel processing to use. Default is "multisession".
#'
#' @return A data frame summarizing all simulation results.
#' Failed simulations are padded with appropriate placeholders to match the structure of successful simulations.
#'
#' @import cmdstanr
#' @import dplyr
#' @import future
#' @import future.apply
#' @import progressr
#'
#' @examples
#' # Define the simulation grid list
#' sim_grid_list <- list(
#'   n = 300,
#'   theta_n = 2,
#'   item_n1 = 6,
#'   item_n2 = 6,
#'   var_ers = c(0.5, 0.3),
#'   var_ars = c(0.5, 0.3),
#'   var_thetas1 = 1,
#'   var_thetas2 = 1,
#'   cor_ers1 = -0.1,
#'   cor_ers2 = 0.1,
#'   cor_thetas = 0.2,
#'   ars_prior = 1,
#'   x_num1 = c(0, 3),
#'   x_num2 = c(0, 3),
#'   chains = 2,
#'   iter = 3000,
#'   warmup = 2000,
#'   adapt_delta = 0.8,
#'   seed = NA,
#'   stan_model = as.character("stan/ARS_ERS.stan"),
#'   init_vals = FALSE,
#'   replication = 1:10
#' )
#'
#' # Create the simulation grid as a data frame
#' sim_grid <- expand.grid(sim_grid_list)
#'
#' # Optional filtering steps
#' sim_grid <- sim_grid[sim_grid$item_n1 == sim_grid$item_n2, ]
#' sim_grid <- sim_grid[sim_grid$x_num1 == sim_grid$x_num2, ]
#' sim_grid <- sim_grid[sim_grid$var_ers == sim_grid$var_ars, ]
#'
#' # Add index column and unique seeds
#' sim_grid$index <- 1:nrow(sim_grid)
#' sim_grid <- sim_grid[, c(ncol(sim_grid), 1:(ncol(sim_grid) - 1))]
#' sim_grid$seed <- sample(1:1e7, size = nrow(sim_grid), replace = FALSE)
#'
#' # Run the simulation function with the grid and save results
#' results_df <- run_simulations_parallel(
#'   sim_grid = sim_grid,
#'   results_path = "results",
#'   parallel_type = "sequential"
#' )
#'
#' # View the summary of simulation results
#' print(results_df)
#'
#' 
#' 
#' @export

sim_fun <- function(sim_grid, results_path, workers = "half", save_all = TRUE, parallel_type = "multisession") {

  # half would mean two cores per replication, quarter would mean four
  if (workers == "half") workers <- future::availableCores() / 2 |> floor()
  if (workers == "quarter") workers <- future::availableCores() / 4 |> floor()

  # Initialize lists to store warnings and errors
  warnings_list <- list()
  errors_list <- list()

  # important if there are multiple models
  # we want to compile them only once
  stan_model_files <- unique(sim_grid$stan_model) |> as.character()

  # Compile all models once in the parent process
  model_files <- unique(sim_grid$stan_model) |> as.character()
  invisible(lapply(model_files, cmdstanr::cmdstan_model))
  
  # write a function to put it into the loop
  get_results <- function(row, results_path) {

    args <- list(
      index      = row["index"] |> as.numeric(),
      n          = row["n"] |> as.numeric(),
      theta_n    = row["theta_n"] |> as.numeric(),
      item_n     = row[grepl("item_n", names(row))] |> unlist() |> as.numeric(),
      var_ers    = row["var_ers"] |> as.numeric(),
      var_ars    = row["var_ars"] |> as.numeric(),
      var_thetas = row[grepl("var_thetas", names(row))] |> unlist() |> as.numeric(),
      cor_ers    = row[grepl("cor_ers", names(row))] |> unlist() |> as.numeric(),
      cor_thetas = row[grepl("cor_thetas", names(row))] |> unlist() |> as.numeric(),
      ars_prior  = row["ars_prior"] |> as.numeric(),
      x_num      = row[grepl("x_num", names(row))] |> unlist() |> as.numeric(),
      chains     = row["chains"] |> as.numeric(),
      iter       = row["iter"] |> as.numeric(),
      warmup     = row["warmup"] |> as.numeric(),
      seed       = row["seed"] |> as.integer(),
      stan_model = row["stan_model"] |> as.character(),
      init_vals  = row["init_vals"] |> as.logical()
    )

    
    # now build a tryCatch to handle warnings and errors
    results <- tryCatch(
      withCallingHandlers({
        do.call(one_simulation, args)
      }, warning = function(w) {
        warnings_list <<- c(conditionMessage(w))
        message(sprintf("Warning during model fitting, index %s: ", row["index"]), 
                conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) {
        errors_list <<- c(conditionMessage(e))
        message(sprintf("Simulation failed, index %s: ", row["index"]), 
                conditionMessage(e))
        return(list(summary = NA, warnings = NULL, error = conditionMessage(e)))
      }
    )
    
    # write warnings
    if (length(warnings_list) > 0) results$warnings <- warnings_list 

    # build a results row with the sim_grid row and the results
    results$summary <- c(row, results$summary) |> as.data.frame() 

    # if there are NAs (non-convergance), be sure it still works
    if (any(is.na(results$summary))) results$summary <- t(results$summary) |> as.data.frame()
    
    # now save this row
    df_row <- results$summary 

    # before continuing, save the results of this replication
    if (save_all == TRUE) save(results, file = sprintf("%s/results_%s.rda", results_path, row["index"]))

    rm(results)

    # garbage collection to save space
    gc()

    # return the important row
    return(df_row)
  }

  # for the progressbar
  handlers("txtprogressbar")

  # sessiontype
  if (parallel_type == "multisession") future::plan(multisession, workers = workers)
  if (parallel_type == "multicore") future::plan(multicore, workers = workers)
  if (parallel_type == "sequential") future::plan(sequential)

  # do the loop with the progressbar and future along the rows of the grid
  with_progress({
    p <- progressor(along = 1:nrow(sim_grid))
    results_df <- future.apply::future_apply(
      sim_grid, 
      MARGIN = 1, 
      FUN = function(row) {
        p(sprintf("Processing a row"))
        get_results(row, results_path = results_path)
      },
      future.seed = TRUE
    )
    Sys.sleep(3)
  })

  # right now in our results_df, we have a list with one row per single replication (so not a df yet)
  # check if there are any failed sims
  # these only have on col more than the sim_grid
  failed_sims <- which(sapply(results_df, function(df) ncol(df) == ncol(sim_grid) + 1))

  # get a simulation row which worked (more cols)
  example_sim_num <- which(sapply(results_df, function(df) ncol(df) > ncol(sim_grid) + 1))[1]
  # get the column names
  example_sim_cols <- colnames(results_df[[example_sim_num]])

  # add NAs to the failed sims based on the working example we extracted
  for (i in failed_sims) {
    tmp <- results_df[[i]]
    missing_cols <- length(example_sim_cols) - ncol(tmp)
    results_df[[i]] <- cbind(tmp, matrix(NA, ncol = missing_cols))
    colnames(results_df[[i]]) <- example_sim_cols
  }

  # add all the rows together for a dataframe
  results_df <- do.call(rbind, results_df)
  rownames(results_df) <- 1:nrow(results_df)

  return(results_df)
}
