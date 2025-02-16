########################################
## run_simulations.R
##
## This script runs the IPCW/AIPCW method across multiple 
## datasets (generated with different random seeds). 
## It provides a wrapper for running experiments efficiently.
########################################

## Load required libraries
suppressMessages(library(BiocManager))   # Bioconductor package manager
suppressMessages(library(BiocParallel))  # Parallel computing support

#' Run a single simulation experiment
#'
#' This function executes the following steps:
#'  1) Generates synthetic survival data based on the specified data-generating process.
#'  2) Applies the IPCW/AIPCW method to estimate lower predictive bounds and empirical coverage.
#'
#' @param seed Random seed for reproducibility.
#' @param data_generating_process A string specifying the data generation method: 
#'   - "exp" for exponential survival data.
#'   - "lognorm" for lognormal survival data.
#'   - "mis" for mis-specified survival data.
#' @param model A string specifying the survival model type: 
#'   - "cox" for Cox proportional hazards model.
#'   - "rsf" for Random Survival Forest (RSF).
#'   - "sl" for Super Learner (SL).
#' @param ... Additional parameters passed to the data generation function.
#' @return A list containing results from the IPCW/AIPCW method, including:
#'   - Empirical coverage metrics.
#'   - Estimated beta indices.
#'   - Predictive lower bounds for test and calibration sets.
#'
simulate_experiment <- function(seed, data_generating_process, model, ...) {
  
  ## Source the helper scripts 
  source("utils_data.R")              # Data generation functions
  source("utils_fit_survival.R")      # Survival model fitting functions
  source("utils_methods_aipcw.R")     # IPCW/AIPCW method functions
  
  
  # Step 1: Generate synthetic survival data based on the specified process
  if (data_generating_process == "exp") {
    sim_data <- simulate_exponential_data(seed = seed, ...)  # Pass additional parameters
  } else if (data_generating_process == "lognorm") {
    sim_data <- simulate_lognorm_data(seed = seed, ...)      # Pass additional parameters
  } else if (data_generating_process == "mis") {
    sim_data <- simulate_mis_data(seed = seed, ...)          # Pass additional parameters
  } else {
    stop("Invalid data generating process. Choose from 'exp', 'lognorm', or 'mis'.")
  }
  
  # Step 2: Apply the IPCW/AIPCW method to compute predictive bounds and coverage
  result <- run_methods_aipcw(
    data = sim_data$data,    # Survival dataset
    T_true = sim_data$T_true, # True event times before censoring
    model = model            # Chosen survival model
  )
  
  return(result)  # Return computed coverage and bounds
}

################# EXAMPLE: RUN MULTIPLE SIMULATIONS ##################
## Run multiple simulation experiments in parallel
##
## This example runs the `simulate_experiment` function 100 times 
## using parallel processing to speed up computation.
##
## Adjust the number of workers (`workers = 12`) based on available CPU cores.
##
## In this example, we generate data using the Exponential process and fit 
## the RSF model.
## We do not pass additional parameters for the data-generating process, 
## the default settings are used.

# Define parallel backend with 12 workers
bbparam <- SnowParam(workers = 12, type = "SOCK")  

# Measure execution time and run simulations
system.time({
  results <- bplapply(1:100, 
                      simulate_experiment, 
                      data_generating_process = 'exp', 
                      model = 'rsf',
                      BPPARAM = bbparam)
})

# Save results to file for later analysis
save(results, file = "result.RData")