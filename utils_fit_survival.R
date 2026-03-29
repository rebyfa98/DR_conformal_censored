########################################
## utils_fit_survival.R
##
## Methods that fit survival and censoring models.
########################################

## Load required libraries
suppressMessages(library(survival))  # Cox models
suppressMessages(library(pec))  # Cox models for prediction
suppressMessages(library(randomForestSRC))  # Random survival forests
suppressMessages(library(SuperLearner))  # Super Learner
suppressMessages(library(survSuperLearner))  # Survival Super Learner

#' Define the default Super Learner library
#'
#' @param include_gam Logical indicating whether `survSL.gam` should be added
#' @return A character vector of Super Learner wrappers
#'
default_sl_library <- function(include_gam = FALSE) {
  base_library <- c(
    "survSL.km",
    "survSL.coxph",
    "survSL.expreg",
    "survSL.weibreg",
    "survSL.loglogreg",
    "survSL.rfsrc"
  )

  if (include_gam) {
    c(base_library, "survSL.gam")
  } else {
    base_library
  }
}

#' Build a prediction time grid restricted to the support of the training data
#'
#' @param observed_times Candidate times from the target sample
#' @param train_times Observed times from the training sample
#' @return A sorted vector of prediction times
#'
make_time_grid <- function(observed_times, train_times) {
  if (length(train_times) == 0) {
    stop("Training set has no observed times for the requested process.")
  }

  grid <- sort(unique(
    observed_times[
      observed_times >= min(train_times) &
        observed_times <= max(train_times)
    ]
  ))

  if (length(grid) == 0) {
    sort(unique(train_times))
  } else {
    grid
  }
}

#' Ensure that predictions are stored as a matrix
#'
#' @param x A vector or matrix of predictions
#' @return A matrix of predictions
#'
ensure_matrix <- function(x) {
  if (is.null(dim(x))) {
    matrix(x, ncol = 1)
  } else {
    x
  }
}

#' Fit survival models based on the specified model type
#'
#' Selects and fits a survival model based on user input.
#'
#' @param model A string specifying the model type: "cox", "rsf", or "sl"
#' @param data_train Training dataset
#' @param data_calib Calibration dataset
#' @param data_test Test dataset
#' @param sl_library Library used for Super Learner
#' @param censor_floor Small positive floor for censoring survival predictions
#' @return Model predictions for the calibration and test datasets
#'
fit_models <- function(
    model,
    data_train,
    data_calib,
    data_test,
    sl_library = default_sl_library(),
    censor_floor = 1e-5
) {
  if (model == "cox") {
    fit_models_cox(
      data_train = data_train,
      data_calib = data_calib,
      data_test = data_test,
      censor_floor = censor_floor
    )
  } else if (model == "rsf") {
    fit_models_rsf(
      data_train = data_train,
      data_calib = data_calib,
      data_test = data_test,
      censor_floor = censor_floor
    )
  } else if (model == "sl") {
    fit_models_sl(
      data_train = data_train,
      data_calib = data_calib,
      data_test = data_test,
      sl_library = sl_library,
      censor_floor = censor_floor
    )
  } else {
    stop("Invalid model type. Choose from 'cox', 'rsf', or 'sl'.")
  }
}

#' Fit Cox proportional hazards survival models
#'
#' Fits separate Cox models for survival and censoring times and
#' predicts survival probabilities at specific time points.
#'
#' @param data_train Training dataset
#' @param data_calib Calibration dataset
#' @param data_test Test dataset
#' @param censor_floor Small positive floor for censoring survival predictions
#' @return A list containing:
#'   \item{pred_time_T}{Prediction time points for T (calibration set)}
#'   \item{pred_surv_T}{Predicted survival probabilities for T (calibration set)}
#'   \item{pred_time_T_test}{Prediction time points for T (test set)}
#'   \item{pred_surv_T_test}{Predicted survival probabilities for T (test set)}
#'   \item{pred_time_C}{Prediction time points for censoring (calibration set)}
#'   \item{pred_surv_C}{Predicted survival probabilities for censoring (calibration set)}
#'
fit_models_cox <- function(
    data_train,
    data_calib,
    data_test,
    censor_floor = 1e-5
) {
  p <- ncol(data_train) - 3  # Number of covariates
  xnames <- colnames(data_train)[seq_len(p)]  # Extract covariate names

  # Observed event and censoring times in the training data
  train_event_times <- sort(unique(data_train$censored_T[data_train$event]))
  train_cens_times <- sort(unique(data_train$censored_T[data_train$event_C]))

  # Fit Cox model for survival times
  fmla_T <- as.formula(
    paste("Surv(censored_T, event) ~", paste(xnames, collapse = " + "))
  )
  mdl_T <- coxph(fmla_T, data = data_train, x = TRUE, method = "breslow")

  # Fit Cox model for censoring times
  fmla_C <- as.formula(
    paste("Surv(censored_T, event_C) ~", paste(xnames, collapse = " + "))
  )
  mdl_C <- coxph(fmla_C, data = data_train, x = TRUE, method = "breslow")

  # Determine prediction time points
  pred_time_T <- make_time_grid(
    observed_times = data_calib$censored_T[data_calib$event],
    train_times = train_event_times
  )
  pred_time_C <- make_time_grid(
    observed_times = data_calib$censored_T[data_calib$event_C],
    train_times = train_cens_times
  )

  pred_surv_T <- ensure_matrix(
    predictSurvProb(mdl_T, newdata = data_calib, times = pred_time_T)
  )
  pred_surv_T_test <- ensure_matrix(
    predictSurvProb(mdl_T, newdata = data_test, times = pred_time_T)
  )
  pred_surv_C <- ensure_matrix(
    predictSurvProb(mdl_C, newdata = data_calib, times = pred_time_C)
  )

  # Return survival probability predictions
  list(
    pred_time_T = pred_time_T,
    pred_surv_T = pred_surv_T,
    pred_time_T_test = pred_time_T,
    pred_surv_T_test = pred_surv_T_test,
    pred_time_C = pred_time_C,
    pred_surv_C = pmax(pred_surv_C, censor_floor)
  )
}

#' Fit Random Survival Forest (RSF) survival models
#'
#' Uses RSF to fit separate models for survival and censoring times and
#' predicts survival probabilities at specific time points.
#'
#' @param data_train Training dataset
#' @param data_calib Calibration dataset
#' @param data_test Test dataset
#' @param censor_floor Small positive floor for censoring survival predictions
#' @return A list containing:
#'   \item{pred_time_T}{Prediction time points for T}
#'   \item{pred_surv_T}{Predicted survival probabilities for T}
#'   \item{pred_time_T_test}{Prediction time points for T (test set)}
#'   \item{pred_surv_T_test}{Predicted survival probabilities for T (test set)}
#'   \item{pred_time_C}{Prediction time points for censoring}
#'   \item{pred_surv_C}{Predicted survival probabilities for censoring}
#'
fit_models_rsf <- function(
    data_train,
    data_calib,
    data_test,
    censor_floor = 1e-5
) {
  p <- ncol(data_train) - 3  # Number of covariates
  xnames <- colnames(data_train)[seq_len(p)]  # Extract covariate names

  # Fit RSF model for survival times
  fmla_T <- as.formula(
    paste("Surv(censored_T, event) ~", paste(xnames, collapse = " + "))
  )
  mdl_T <- rfsrc(fmla_T, data = data_train, ntime = NULL)

  # Fit RSF model for censoring times
  fmla_C <- as.formula(
    paste("Surv(censored_T, event_C) ~", paste(xnames, collapse = " + "))
  )
  mdl_C <- rfsrc(fmla_C, data = data_train, ntime = NULL)

  # Generate predictions for calibration and test sets
  pred_T <- predict(mdl_T, data_calib)
  pred_T_test <- predict(mdl_T, data_test)
  pred_C <- predict(mdl_C, data_calib)

  train_event_times <- sort(unique(data_train$censored_T[data_train$event]))
  train_cens_times <- sort(unique(data_train$censored_T[data_train$event_C]))

  # Align prediction grids with the current AIPCW implementation
  pred_time_T <- if (length(train_event_times) == ncol(pred_T$survival)) {
    train_event_times
  } else {
    pred_T$time.interest
  }
  pred_time_C <- if (length(train_cens_times) == ncol(pred_C$survival)) {
    train_cens_times
  } else {
    pred_C$time.interest
  }

  # Return survival probability predictions
  list(
    pred_time_T = pred_time_T,
    pred_surv_T = ensure_matrix(pred_T$survival),
    pred_time_T_test = pred_time_T,
    pred_surv_T_test = ensure_matrix(pred_T_test$survival),
    pred_time_C = pred_time_C,
    pred_surv_C = pmax(ensure_matrix(pred_C$survival), censor_floor)
  )
}

#' Fit Super Learner survival models
#'
#' Uses Super Learner (SL) to fit an ensemble model for survival prediction.
#'
#' @param data_train Training dataset
#' @param data_calib Calibration dataset
#' @param data_test Test dataset
#' @param sl_library Library used for Super Learner
#' @param censor_floor Small positive floor for censoring survival predictions
#' @return A list containing:
#'   \item{pred_time_T}{Prediction time points for T}
#'   \item{pred_surv_T}{Predicted survival probabilities for T}
#'   \item{pred_time_T_test}{Prediction time points for T (test set)}
#'   \item{pred_surv_T_test}{Predicted survival probabilities for T (test set)}
#'   \item{pred_time_C}{Prediction time points for censoring}
#'   \item{pred_surv_C}{Predicted survival probabilities for censoring}
#'
fit_models_sl <- function(
    data_train,
    data_calib,
    data_test,
    sl_library = default_sl_library(),
    censor_floor = 1e-5
) {
  p <- ncol(data_train) - 3  # Number of covariates

  # Observed event and censoring times in the training data
  train_event_times <- sort(unique(data_train$censored_T[data_train$event]))
  train_cens_times <- sort(unique(data_train$censored_T[data_train$event_C]))

  # Determine prediction time points
  pred_time_T <- make_time_grid(
    observed_times = data_calib$censored_T[data_calib$event],
    train_times = train_event_times
  )
  pred_time_T_test <- make_time_grid(
    observed_times = data_test$censored_T[data_test$event],
    train_times = train_event_times
  )
  pred_time_C <- make_time_grid(
    observed_times = data_calib$censored_T[data_calib$event_C],
    train_times = train_cens_times
  )

  # Fit Super Learner model for survival
  fit <- survSuperLearner(
    time = data_train$censored_T,
    event = data_train$event,
    X = data_train[, seq_len(p), drop = FALSE],
    newX = data_calib[, seq_len(p), drop = FALSE],
    new.times = pred_time_T,
    event.SL.library = sl_library,
    cens.SL.library = sl_library,
    verbose = FALSE
  )

  pred_T_test <- predict.survSuperLearner(
    fit,
    newdata = data_test[, seq_len(p), drop = FALSE],
    new.times = pred_time_T_test,
    onlySL = TRUE
  )
  pred_C <- predict.survSuperLearner(
    fit,
    newdata = data_calib[, seq_len(p), drop = FALSE],
    new.times = pred_time_C,
    onlySL = TRUE
  )

  # Return survival probability predictions
  list(
    pred_time_T = pred_time_T,
    pred_surv_T = ensure_matrix(fit$event.SL.predict),
    pred_time_T_test = pred_time_T_test,
    pred_surv_T_test = ensure_matrix(pred_T_test$event.SL.predict),
    pred_time_C = pred_time_C,
    pred_surv_C = pmax(ensure_matrix(pred_C$cens.SL.predict), censor_floor)
  )
}
