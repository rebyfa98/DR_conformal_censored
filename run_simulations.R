########################################
## run_simulations.R
##
## This script runs the IPCW/AIPCW method across multiple
## datasets (generated with different random seeds).
## It provides a wrapper for running experiments efficiently.
########################################

## Load required libraries
suppressMessages(library(BiocParallel))  # Parallel computing support

#' Run a single simulation experiment
#'
#' This function executes the following steps:
#'  1) Generates synthetic survival data based on the specified setting.
#'  2) Applies the IPCW/AIPCW method to estimate lower predictive bounds and empirical coverage.
#'
#' @param seed Random seed for reproducibility.
#' @param setting A string specifying the synthetic setting:
#'   - "setting1", ..., "setting6".
#' @param model A string specifying the survival model type:
#'   - "cox" for Cox proportional hazards model.
#'   - "rsf" for Random Survival Forest (RSF).
#'   - "sl" for Super Learner (SL).
#' @param n_train Number of training samples.
#' @param n_calib Number of calibration samples.
#' @param n_test Number of test samples.
#' @param alpha Miscoverage level.
#' @param beta_grid Quantile grid used by the conformal procedure.
#' @param sl_library Super Learner library used when `model = "sl"`.
#' @param censor_floor Lower floor for estimated censoring survival probabilities.
#' @return A list containing results from the IPCW/AIPCW method, including:
#'   - Empirical coverage metrics.
#'   - Estimated beta indices.
#'   - Predictive lower bounds for test and calibration sets.
#'
simulate_experiment <- function(
    seed,
    setting,
    model,
    n_train = 1000,
    n_calib = 1000,
    n_test = 1000,
    alpha = 0.1,
    beta_grid = seq(0, 1, by = 0.001),
    sl_library = default_sl_library(),
    censor_floor = 1e-5
) {

  ## Source the helper scripts
  source("utils_data.R")              # Data generation functions
  source("utils_fit_survival.R")      # Survival model fitting functions
  source("utils_methods_aipcw.R")     # IPCW/AIPCW method functions

  # Step 1: Generate synthetic survival data based on the specified setting
  sim_data <- simulate_data(
    setting = setting,
    seed = seed,
    n = n_train + n_calib + n_test
  )

  # Step 2: Apply the IPCW/AIPCW method to compute predictive bounds and coverage
  result <- run_methods_aipcw(
    data = sim_data$data,      # Survival dataset
    T_true = sim_data$T_true,  # True event times before censoring
    model = model,             # Chosen survival model
    n_train = n_train,
    n_calib = n_calib,
    alpha = alpha,
    beta_grid = beta_grid,
    sl_library = sl_library,
    censor_floor = censor_floor
  )

  return(result)  # Return computed coverage and bounds
}

################# EXAMPLE: RUN MULTIPLE SIMULATIONS ##################
## Run multiple simulation experiments in parallel
##
## This example runs the `simulate_experiment` function 100 times
## using parallel processing to speed up computation.
##
## Adjust the number of workers (`workers = 20`) based on available CPU cores.
##
## In this example, we generate data using setting 1
## and fit the RSF model.

# Define parallel backend with 20 workers
bbparam <- SnowParam(workers = 20, type = "SOCK")

# Measure execution time and run simulations
system.time({
  results <- bplapply(
    1:100,
    simulate_experiment,
    setting = "setting1",
    model = "sl",
    BPPARAM = bbparam
  )
})

# Save results to file for later analysis
save(results, file = "result.RData")