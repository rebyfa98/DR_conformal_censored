########################################
## utils_methods_aipcw.R
##
## Methods that run the procedure on the data
## and output coverage, lower predictive bounds, etc.
########################################

#' Compute the estimated survival probability at a given time
#'
#' Retrieves the estimated survival probability for an individual `i` at time `t`,
#' based on the predicted survival curve.
#'
#' @param i Index of the individual in the calibration dataset
#' @param t Time at which survival probability is evaluated
#' @param pred_time Vector of time points for which survival probabilities were estimated
#' @param pred_surv Matrix of survival probabilities (rows correspond to individuals, columns to time points)
#' @return Estimated survival probability at time `t` for individual `i`
#'
surv_curve <- function(
    i, t, pred_time, pred_surv
) {
  # If the requested time t is before the first prediction time, return survival probability 1
  if (t < min(pred_time)) {
    return(1)
  } else {
    # Find the largest time point in pred_time that is less than or equal to t
    index <- max(which(pred_time <= t))
    return(pred_surv[i, index])  # Return the corresponding survival probability
  }
}

#' Compute the estimated survival quantile for a given probability level
#'
#' Finds the smallest time `t` such that the survival probability is at most `1 - beta_grid[j]`
#' for an individual `i`. This corresponds to the estimated quantile of the survival distribution.
#'
#' @param i Index of the individual in the calibration dataset
#' @param j Index in the beta_grid, corresponding to a probability threshold
#' @param pred_time Vector of time points for which survival probabilities were estimated
#' @param pred_surv Matrix of survival probabilities (rows correspond to individuals, columns to time points)
#' @param beta_grid Vector of probability thresholds for quantile computation
#' @return Estimated survival quantile for individual `i` at probability level `beta_grid[j]`
#'
surv_quantile <- function(
    i, j, pred_time, pred_surv, beta_grid
) {
  # If the probability threshold is greater than the final survival probability, return the max time
  if (beta_grid[j] > 1 - pred_surv[i, length(pred_time)]) {
    return(max(pred_time))
  } else {
    # Find the smallest time index where survival probability is at most 1 - beta_grid[j]
    return(pred_time[min(which(1 - pred_surv[i, ] >= beta_grid[j]))])
  }
}

#' Run IPCW and AIPCW methods to compute lower predictive bounds and coverage
#'
#' This function implements Inverse Probability of Censoring Weighting (IPCW) and Augmented IPCW (AIPCW) 
#' to estimate predictive lower bounds for survival times and empirical coverage.
#'
#' @param data A data frame containing survival and censoring information.
#' @param T_true Vector of true event times corresponding to `data`.
#' @param model A string specifying the survival model type ("cox", "rsf", or "sl").
#' @param n_train Number of training samples.
#' @param n_calib Number of calibration samples.
#' @param alpha Desired coverage level (e.g., 0.1 means 90% coverage).
#' @param beta_grid A numeric sequence defining quantile grid points for prediction.
#' @return A list containing empirical coverage rates, estimated beta indices, and lower bounds.
#'
run_methods_aipcw <- function(
    data,
    T_true, 
    model,
    n_train = 1000,
    n_calib = 1000,
    alpha = 0.1, 
    beta_grid = seq(0, 1, by = 0.001)
) {
  
  n <- nrow(data)
  n_test <- n - n_train - n_calib  # Compute number of test samples
  
  # Split data into training, calibration, and test sets
  data_train <- data[1:n_train, ]
  data_calib <- data[(n_train + 1):(n_train + n_calib), ]
  data_test  <- data[(n_train + n_calib + 1):n, ]
  T_calib <- T_true[(n_train + 1):(n_train + n_calib)]
  T_test  <- T_true[(n_train + n_calib + 1):n] 
  
  # Fit the specified survival model and obtain predictions
  pred <- fit_models(model, data_train, data_calib, data_test)
  pred_time_T <- pred$pred_time_T
  pred_surv_T <- pred$pred_surv_T
  pred_time_T_test <- pred$pred_time_T_test
  pred_surv_T_test <- pred$pred_surv_T_test 
  pred_time_C <- pred$pred_time_C
  pred_surv_C <- pred$pred_surv_C
  
  # Compute matrix M_C, used for the censoring process
  u <- sort(unique(data_calib$censored_T[data_calib$event_C]))  # Unique censoring times
  K <- length(u)
  M_C <- matrix(nrow = n_calib, ncol = K + 1)
  
  for (i in 1:n_calib) {
    M_C[i, 1] <- (data_calib$censored_T[i] <= 0 & data_calib$event_C[i]) + 
      log(surv_curve(i, min(0, data_calib$censored_T[i]), pred_time_C, pred_surv_C))
    for (k in 2:(K + 1)) {
      M_C[i, k] <- (data_calib$censored_T[i] <= u[k - 1] & data_calib$event_C[i]) + 
        log(surv_curve(i, min(u[k - 1], data_calib$censored_T[i]), pred_time_C, pred_surv_C))
    }
  }
  
  # Compute matrix D_C based on M_C
  D_C <- matrix(nrow = n_calib, ncol = K)
  for (i in 1:n_calib) {
    for (k in 2:(K + 1)) {
      D_C[i, k - 1] <- M_C[i, k] - M_C[i, k - 1]
    }
  }
  
  # Compute hat_W and hat_Pi
  hat_W <- 1
  hat_Pi <- 1
  j <- 1
  while ((hat_W[j] >= 0 | hat_W[j] + hat_Pi[j] >= 0) & j<=length(beta_grid)) {
    
    ## hat_W
    num <- numeric(n_calib)
    for (i in 1:n_calib) {
      if (!data_calib$event[i]) {
        num[i] <- 0
      } else {
        num[i] <- ((data_calib$censored_T[i] >= 
                      surv_quantile(i,j, pred_time_T, pred_surv_T,beta_grid)) -
                     (1 - alpha)) /
          surv_curve(i, data_calib$censored_T[i], pred_time_C, pred_surv_C)
      }
    }
    hat_W <- c(hat_W, mean(num))
    
    ## hat_Pi
    integral <- numeric(n_calib)
    for (i in 1:n_calib) {
      int_i <- numeric(K)
      for (k in 1:K) {
        if (D_C[i,k] != 0) {
          hat_eta <- surv_curve(i, 
                                max(surv_quantile(i,j, pred_time_T, 
                                                  pred_surv_T, beta_grid), u[k]), 
                                pred_time_T, pred_surv_T) / 
            surv_curve(i, u[k],pred_time_T, pred_surv_T)
          int_i[k] <- (hat_eta - (1 - alpha)) /
            surv_curve(i, u[k], pred_time_C, pred_surv_C) * D_C[i, k]
          if (is.na(int_i[k])) int_i[k] <- 0
        } else {
          int_i[k] <- 0
        }
      }
      integral[i] <- sum(int_i)
    }
    hat_Pi <- c(hat_Pi, mean(integral))
    j <- j + 1
  }
  
  # Remove the first placeholder
  hat_W  <- hat_W[-1]
  hat_Pi <- hat_Pi[-1]
  
  # Select beta index for IPCW method
  j <- 1
  while (hat_W[j] >= 0) {
    j <- j + 1
  }
  hat_beta_index_ipcw <- max(j - 1, 1) # Ensure at least index 1
  
  # Compute IPCW lower bounds
  lower_bnd_calib_ipcw <- sapply(1:n_calib, surv_quantile, 
                                 hat_beta_index_ipcw, pred_time_T, 
                                 pred_surv_T, beta_grid)
  lower_bnd_test_ipcw <- sapply(1:n_test, surv_quantile, 
                                hat_beta_index_ipcw, pred_time_T_test, 
                                pred_surv_T_test, beta_grid)
  
  # Select beta index for AIPCW method
  j <- 1
  while ((hat_W + hat_Pi)[j] >= 0) {
    j <- j + 1
  }
  hat_beta_index_aipcw <- max(j - 1, 1)  # Ensure at least index 1
  
  # Compute AIPCW lower bounds
  lower_bnd_calib_aipcw  <- sapply(1:n_calib, surv_quantile, 
                                   hat_beta_index_aipcw, pred_time_T, 
                                   pred_surv_T, beta_grid)
  lower_bnd_test_aipcw  <- sapply(1:n_test, surv_quantile, 
                                  hat_beta_index_aipcw, pred_time_T_test, 
                                  pred_surv_T_test, beta_grid)
  
  # Compute empirical coverage for calibration and test sets
  list(
    emp_coverage_calib_ipcw = mean(T_calib >= lower_bnd_calib_ipcw),
    emp_coverage_calib_aipcw = mean(T_calib >= lower_bnd_calib_aipcw),
    emp_coverage_test_ipcw = mean(T_test >= lower_bnd_test_ipcw),
    emp_coverage_test_aipcw = mean(T_test >= lower_bnd_test_aipcw),
    hat_beta_index_ipcw = hat_beta_index_ipcw,
    hat_beta_index_aipcw = hat_beta_index_aipcw,
    lower_bnd_test_ipcw = lower_bnd_test_ipcw,
    lower_bnd_test_aipcw = lower_bnd_test_aipcw
  )
}

#' Run OR and COR methods to compute lower predictive bounds and coverage
#'
#' This function implements Outcome Regression (OR) and Calibrated Outcome
#' Regression (COR) to estimate predictive lower bounds for survival times 
#' and empirical coverage.
#'
#' @param data A data frame containing survival and censoring information.
#' @param T_true Vector of true event times corresponding to `data`.
#' @param model A string specifying the survival model type ("cox", "rsf", or "sl").
#' @param n_train Number of training samples.
#' @param n_calib Number of calibration samples.
#' @param alpha Desired coverage level (e.g., 0.1 means 90% coverage).
#' @param beta_grid A numeric sequence defining quantile grid points for prediction.
#' @return A list containing empirical coverage rates, estimated beta indices, and lower bounds.
#'
run_methods_cor <- function(
    data,
    T_true, 
    model,
    n_train = 1000,
    n_calib = 1000,
    alpha = 0.1, 
    beta_grid = seq(0, 1, by = 0.001)
) {
  
  n <- nrow(data)
  n_test <- n - n_train - n_calib  # Compute number of test samples
  
  # Split data into training, calibration, and test sets
  data_train <- data[1:n_train, ]
  data_calib <- data[(n_train + 1):(n_train + n_calib), ]
  data_test  <- data[(n_train + n_calib + 1):n, ]
  T_calib <- T_true[(n_train + 1):(n_train + n_calib)]
  T_test  <- T_true[(n_train + n_calib + 1):n] 
  
  # Fit the specified survival model and obtain predictions
  pred <- fit_models(model, data_train, data_calib, data_test)
  pred_time_T <- pred$pred_time_T
  pred_surv_T <- pred$pred_surv_T
  pred_time_T_test <- pred$pred_time_T_test
  pred_surv_T_test <- pred$pred_surv_T_test 
  pred_time_C <- pred$pred_time_C
  pred_surv_C <- pred$pred_surv_C
  
  # OR method: the estimated LPB is the predicted quantile at level alpha
  alpha_index <- which.min(abs(beta_grid - alpha))
  
  # COR method
  hat_P <- 1
  j <- 1
  while (hat_P[j] >= 1-alpha){
    num <- numeric(n_calib)
    for (i in 1:n_calib) {
      num[i] <- surv_curve(i, surv_quantile(i, j, pred_time_T, pred_surv_T,beta_grid),
                           pred_time_T, pred_surv_T)
    }
    hat_P <- c(hat_P, mean(num))
    j <- j + 1 
  }
  
  # Remove the first placeholder
  hat_P <- hat_P[-1]
  
  # Select beta index
  hat_beta_index <- max(j - 1, 1) # Ensure at least index 1
  hat_beta <- beta_grid[hat_beta_index]
  
  # Compute OR and COR lower bounds
  lower_bnd_calib_or  <- sapply(1:n_calib, surv_quantile, 
                                alpha_index, pred_time_T, 
                                pred_surv_T, beta_grid)
  lower_bnd_test_or  <- sapply(1:n_test, surv_quantile, 
                               alpha_index, pred_time_T_test, 
                               pred_surv_T_test, beta_grid)
  lower_bnd_calib_cor  <- sapply(1:n_calib, surv_quantile, 
                                 hat_beta_index, pred_time_T, 
                                 pred_surv_T, beta_grid)
  lower_bnd_test_cor  <- sapply(1:n_test, surv_quantile, 
                                hat_beta_index, pred_time_T_test, 
                                pred_surv_T_test, beta_grid)
  
  # Compute empirical coverage for calibration and test sets
  list(
    emp_coverage_calib_or = mean(T_calib >= lower_bnd_calib_or),
    emp_coverage_calib_cor = mean(T_calib >= lower_bnd_calib_cor),
    emp_coverage_test_or = mean(T_test >= lower_bnd_test_or),
    emp_coverage_test_cor = mean(T_test >= lower_bnd_test_cor),
    hat_beta_index = hat_beta_index,
    lower_bnd_test_or = lower_bnd_test_or,
    lower_bnd_test_cor = lower_bnd_test_cor
  )
}


