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
#' @param i Index of the individual in the dataset
#' @param t Time at which survival probability is evaluated
#' @param pred_time Vector of time points for which survival probabilities were estimated
#' @param pred_surv Matrix of survival probabilities
#' @return Estimated survival probability at time `t` for individual `i`
#'
surv_curve <- function(i, t, pred_time, pred_surv) {
  if (length(pred_time) == 0) {
    stop("Prediction time grid is empty.")
  }

  # If the requested time t is before the first prediction time, return survival probability 1
  if (t < min(pred_time)) {
    1
  } else {
    # Find the largest time point in pred_time that is less than or equal to t
    index <- max(which(pred_time <= t))
    pred_surv[i, index]
  }
}

#' Compute the estimated survival quantile for a given probability level
#'
#' Finds the smallest time `t` such that the survival probability is at most `1 - beta_grid[j]`
#' for an individual `i`.
#'
#' @param i Index of the individual in the dataset
#' @param j Index in the beta_grid
#' @param pred_time Vector of time points for which survival probabilities were estimated
#' @param pred_surv Matrix of survival probabilities
#' @param beta_grid Vector of probability thresholds for quantile computation
#' @return Estimated survival quantile for individual `i` at probability level `beta_grid[j]`
#'
surv_quantile <- function(i, j, pred_time, pred_surv, beta_grid) {
  if (length(pred_time) == 0) {
    stop("Prediction time grid is empty.")
  }

  # If the probability threshold is greater than the final survival probability, return the max time
  if (beta_grid[j] > 1 - pred_surv[i, length(pred_time)]) {
    max(pred_time)
  } else {
    # Find the smallest time index where survival probability is at most 1 - beta_grid[j]
    crossing <- which(1 - pred_surv[i, ] >= beta_grid[j])
    if (length(crossing) == 0) {
      max(pred_time)
    } else {
      pred_time[min(crossing)]
    }
  }
}

#' Compute the estimated IPCW and AIPCW score terms over the beta grid
#'
#' @param data_calib Calibration dataset
#' @param pred_time_T Prediction times for the event process
#' @param pred_surv_T Predicted survival probabilities for the event process
#' @param pred_time_C Prediction times for the censoring process
#' @param pred_surv_C Predicted survival probabilities for the censoring process
#' @param alpha Desired miscoverage level
#' @param beta_grid Quantile grid
#' @return A list containing the estimated IPW and influence-function terms
#'
compute_hat_terms <- function(
    data_calib,
    pred_time_T,
    pred_surv_T,
    pred_time_C,
    pred_surv_C,
    alpha,
    beta_grid
) {
  # Unique censoring times in the calibration data
  u <- sort(unique(data_calib$censored_T[data_calib$event_C]))
  if (length(u) == 0) {
    stop("Calibration set has no censoring times.")
  }

  n_calib <- nrow(data_calib)
  K <- length(u)

  # Compute matrix M_C, used for the censoring process
  M_C <- matrix(0, nrow = n_calib, ncol = K + 1)
  for (i in seq_len(n_calib)) {
    M_C[i, 1] <- (data_calib$censored_T[i] <= 0 && data_calib$event_C[i]) +
      log(surv_curve(i, min(0, data_calib$censored_T[i]), pred_time_C, pred_surv_C))
    for (k in 2:(K + 1)) {
      M_C[i, k] <- (data_calib$censored_T[i] <= u[k - 1] && data_calib$event_C[i]) +
        log(surv_curve(i, min(u[k - 1], data_calib$censored_T[i]), pred_time_C, pred_surv_C))
    }
  }

  # Compute matrix D_C based on M_C
  D_C <- M_C[, -1, drop = FALSE] - M_C[, -ncol(M_C), drop = FALSE]

  hat_P <- numeric(0)
  hat_IF <- numeric(0)

  for (j in seq_along(beta_grid)) {
    ## hat_P
    num <- numeric(n_calib)
    for (i in seq_len(n_calib)) {
      if (!data_calib$event[i]) {
        num[i] <- 0
      } else {
        q_ij <- surv_quantile(i, j, pred_time_T, pred_surv_T, beta_grid)
        num[i] <- ((data_calib$censored_T[i] >= q_ij) - (1 - alpha)) /
          surv_curve(i, data_calib$censored_T[i], pred_time_C, pred_surv_C)
      }
    }
    hat_P[j] <- mean(num)

    ## hat_IF
    integral <- numeric(n_calib)
    for (i in seq_len(n_calib)) {
      int_i <- numeric(K)
      q_ij <- surv_quantile(i, j, pred_time_T, pred_surv_T, beta_grid)
      for (k in seq_len(K)) {
        if (D_C[i, k] != 0) {
          hat_Q <- surv_curve(
            i,
            max(q_ij, u[k]),
            pred_time_T,
            pred_surv_T
          ) / surv_curve(i, u[k], pred_time_T, pred_surv_T)

          int_i[k] <- (hat_Q - (1 - alpha)) /
            surv_curve(i, u[k], pred_time_C, pred_surv_C) *
            D_C[i, k]

          if (is.na(int_i[k])) {
            int_i[k] <- 0
          }
        }
      }
      integral[i] <- sum(int_i)
    }
    hat_IF[j] <- mean(integral)

    if (hat_P[j] < 0 && (hat_P[j] + hat_IF[j]) < 0) {
      break
    }
  }

  list(
    hat_P = hat_P,
    hat_IF = hat_IF
  )
}

#' Find the last nonnegative entry of a numeric vector
#'
#' @param x Numeric vector
#' @return The last index for which `x` is nonnegative, or 1 if none exists
#'
last_nonnegative_index <- function(x) {
  idx <- which(x >= 0)
  if (length(idx) == 0) {
    1L
  } else {
    max(idx)
  }
}

#' Run IPCW and AIPCW methods to compute lower predictive bounds and coverage
#'
#' This function implements Inverse Probability of Censoring Weighting (IPCW) and
#' Augmented IPCW (AIPCW) to estimate predictive lower bounds for survival times
#' and empirical coverage.
#'
#' @param data A data frame containing survival and censoring information
#' @param T_true Vector of true event times corresponding to `data`
#' @param model A string specifying the survival model type ("cox", "rsf", or "sl")
#' @param n_train Number of training samples
#' @param n_calib Number of calibration samples
#' @param alpha Desired coverage level
#' @param beta_grid A numeric sequence defining quantile grid points for prediction
#' @param sl_library Library used for Super Learner
#' @param censor_floor Small positive floor for censoring survival predictions
#' @return A list containing empirical coverage rates, estimated beta indices, and lower bounds
#'
run_methods_aipcw <- function(
    data,
    T_true,
    model,
    n_train = 1000,
    n_calib = 1000,
    alpha = 0.1,
    beta_grid = seq(0, 1, by = 0.001),
    sl_library = default_sl_library(),
    censor_floor = 1e-5
) {
  n <- nrow(data)
  n_test <- n - n_train - n_calib  # Compute number of test samples

  # Split data into training, calibration, and test sets
  data_train <- data[seq_len(n_train), , drop = FALSE]
  data_calib <- data[(n_train + 1):(n_train + n_calib), , drop = FALSE]
  data_test <- data[(n_train + n_calib + 1):n, , drop = FALSE]

  T_calib <- T_true[(n_train + 1):(n_train + n_calib)]
  T_test <- T_true[(n_train + n_calib + 1):n]

  # Fit the specified survival model and obtain predictions
  pred <- fit_models(
    model = model,
    data_train = data_train,
    data_calib = data_calib,
    data_test = data_test,
    sl_library = sl_library,
    censor_floor = censor_floor
  )

  # Compute estimated IPW and influence-function terms
  hat_terms <- compute_hat_terms(
    data_calib = data_calib,
    pred_time_T = pred$pred_time_T,
    pred_surv_T = pred$pred_surv_T,
    pred_time_C = pred$pred_time_C,
    pred_surv_C = pred$pred_surv_C,
    alpha = alpha,
    beta_grid = beta_grid
  )

  # Select beta indices for IPCW and AIPCW
  hat_beta_index_ipcw <- last_nonnegative_index(hat_terms$hat_P)
  hat_beta_index_aipcw <- last_nonnegative_index(hat_terms$hat_P + hat_terms$hat_IF)

  # Compute lower bounds
  lower_bnd_calib_ipcw <- sapply(
    seq_len(n_calib),
    surv_quantile,
    j = hat_beta_index_ipcw,
    pred_time = pred$pred_time_T,
    pred_surv = pred$pred_surv_T,
    beta_grid = beta_grid
  )
  lower_bnd_test_ipcw <- sapply(
    seq_len(n_test),
    surv_quantile,
    j = hat_beta_index_ipcw,
    pred_time = pred$pred_time_T_test,
    pred_surv = pred$pred_surv_T_test,
    beta_grid = beta_grid
  )

  lower_bnd_calib_aipcw <- sapply(
    seq_len(n_calib),
    surv_quantile,
    j = hat_beta_index_aipcw,
    pred_time = pred$pred_time_T,
    pred_surv = pred$pred_surv_T,
    beta_grid = beta_grid
  )
  lower_bnd_test_aipcw <- sapply(
    seq_len(n_test),
    surv_quantile,
    j = hat_beta_index_aipcw,
    pred_time = pred$pred_time_T_test,
    pred_surv = pred$pred_surv_T_test,
    beta_grid = beta_grid
  )

  # Compute empirical coverage for calibration and test sets
  list(
    emp_coverage_calib_ipcw = mean(T_calib >= lower_bnd_calib_ipcw),
    emp_coverage_calib_aipcw = mean(T_calib >= lower_bnd_calib_aipcw),
    emp_coverage_test_ipcw = mean(T_test >= lower_bnd_test_ipcw),
    emp_coverage_test_aipcw = mean(T_test >= lower_bnd_test_aipcw),
    hat_beta_index_ipcw = hat_beta_index_ipcw,
    hat_beta_index_aipcw = hat_beta_index_aipcw,
    lower_bnd_test_ipcw = lower_bnd_test_ipcw,
    lower_bnd_test_aipcw = lower_bnd_test_aipcw,
    hat_P = hat_terms$hat_P,
    hat_IF = hat_terms$hat_IF
  )
}

#' Run OR and COR methods to compute lower predictive bounds and coverage
#'
#' This function implements Outcome Regression (OR) and Calibrated Outcome
#' Regression (COR) to estimate predictive lower bounds for survival times
#' and empirical coverage.
#'
#' @param data A data frame containing survival and censoring information
#' @param T_true Vector of true event times corresponding to `data`
#' @param model A string specifying the survival model type ("cox", "rsf", or "sl")
#' @param n_train Number of training samples
#' @param n_calib Number of calibration samples
#' @param alpha Desired coverage level
#' @param beta_grid A numeric sequence defining quantile grid points for prediction
#' @param sl_library Library used for Super Learner
#' @param censor_floor Small positive floor for censoring survival predictions
#' @return A list containing empirical coverage rates, estimated beta indices, and lower bounds
#'
run_methods_cor <- function(
    data,
    T_true,
    model,
    n_train = 1000,
    n_calib = 1000,
    alpha = 0.1,
    beta_grid = seq(0, 1, by = 0.001),
    sl_library = default_sl_library(),
    censor_floor = 1e-5
) {
  n <- nrow(data)
  n_test <- n - n_train - n_calib  # Compute number of test samples

  # Split data into training, calibration, and test sets
  data_train <- data[seq_len(n_train), , drop = FALSE]
  data_calib <- data[(n_train + 1):(n_train + n_calib), , drop = FALSE]
  data_test <- data[(n_train + n_calib + 1):n, , drop = FALSE]

  T_calib <- T_true[(n_train + 1):(n_train + n_calib)]
  T_test <- T_true[(n_train + n_calib + 1):n]

  # Fit the specified survival model and obtain predictions
  pred <- fit_models(
    model = model,
    data_train = data_train,
    data_calib = data_calib,
    data_test = data_test,
    sl_library = sl_library,
    censor_floor = censor_floor
  )

  # OR method: the estimated LPB is the predicted quantile at level alpha
  alpha_index <- which.min(abs(beta_grid - alpha))

  # COR method
  hat_P <- numeric(0)
  for (j in seq_along(beta_grid)) {
    num <- numeric(n_calib)
    for (i in seq_len(n_calib)) {
      q_ij <- surv_quantile(i, j, pred$pred_time_T, pred$pred_surv_T, beta_grid)
      num[i] <- surv_curve(i, q_ij, pred$pred_time_T, pred$pred_surv_T)
    }
    hat_P[j] <- mean(num)

    if (hat_P[j] < (1 - alpha)) {
      break
    }
  }

  hat_beta_index <- last_nonnegative_index(hat_P - (1 - alpha))

  # Compute OR and COR lower bounds
  lower_bnd_calib_or <- sapply(
    seq_len(n_calib),
    surv_quantile,
    j = alpha_index,
    pred_time = pred$pred_time_T,
    pred_surv = pred$pred_surv_T,
    beta_grid = beta_grid
  )
  lower_bnd_test_or <- sapply(
    seq_len(n_test),
    surv_quantile,
    j = alpha_index,
    pred_time = pred$pred_time_T_test,
    pred_surv = pred$pred_surv_T_test,
    beta_grid = beta_grid
  )
  lower_bnd_calib_cor <- sapply(
    seq_len(n_calib),
    surv_quantile,
    j = hat_beta_index,
    pred_time = pred$pred_time_T,
    pred_surv = pred$pred_surv_T,
    beta_grid = beta_grid
  )
  lower_bnd_test_cor <- sapply(
    seq_len(n_test),
    surv_quantile,
    j = hat_beta_index,
    pred_time = pred$pred_time_T_test,
    pred_surv = pred$pred_surv_T_test,
    beta_grid = beta_grid
  )

  # Compute empirical coverage for calibration and test sets
  list(
    emp_coverage_calib_or = mean(T_calib >= lower_bnd_calib_or),
    emp_coverage_calib_cor = mean(T_calib >= lower_bnd_calib_cor),
    emp_coverage_test_or = mean(T_test >= lower_bnd_test_or),
    emp_coverage_test_cor = mean(T_test >= lower_bnd_test_cor),
    hat_beta_index = hat_beta_index,
    lower_bnd_test_or = lower_bnd_test_or,
    lower_bnd_test_cor = lower_bnd_test_cor,
    hat_P = hat_P
  )
}
