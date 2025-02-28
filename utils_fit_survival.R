########################################
## utils_method.R
##
## Methods that fit survival and censoring models.
########################################

## Load required libraries
suppressMessages(library(survival))  # Cox models
suppressMessages(library(pec)) # Cox models for prediction
suppressMessages(library(randomForestSRC))  # Random survival forests
suppressMessages(library(SuperLearner)) # Super Learner
suppressMessages(library(survSuperLearner)) # Super Learner

#' Fit survival models based on the specified model type
#'
#' Selects and fits a survival model based on user input.
#'
#' @param model A string specifying the model type: "cox", "rsf", or "sl"
#' @param data_train Training dataset
#' @param data_calib Calibration dataset
#' @param data_test Test dataset
#' @return Model predictions for the test dataset
#'
fit_models <- function(model, data_train, data_calib, data_test) {
  if (model == "cox") {
    pred <- fit_models_cox(data_train, data_calib, data_test)
  } else if (model == "rsf") {
    pred <- fit_models_rfsrc(data_train, data_calib, data_test)
  } else if (model == "sl") {
    pred <- fit_models_sl(data_train, data_calib, data_test)
  } else {
    stop("Invalid model type. Choose from 'cox', 'rsf', or 'sl'.")
  }
  
  return(pred)
}

#' Fit Cox proportional hazards survival models
#'
#' Fits separate Cox models for survival and censoring times and 
#' predicts survival probabilities at specific time points.
#'
#' @param data_train Training dataset
#' @param data_calib Calibration dataset
#' @param data_test Test dataset
#' @return A list containing:
#'   \item{pred_time_T}{Prediction time points for T (calibration set)}
#'   \item{pred_surv_T}{Predicted survival probabilities for T (calibration set)}
#'   \item{pred_time_T_test}{Prediction time points for T (test set)}
#'   \item{pred_surv_T_test}{Predicted survival probabilities for T (test set)}
#'   \item{pred_time_C}{Prediction time points for censoring (calibration set)}
#'   \item{pred_surv_C}{Predicted survival probabilities for censoring (calibration set)}
#'
fit_models_cox <- function(
    data_train, data_calib, data_test
) {
  p <- ncol(data_train) - 3  # Number of covariates
  xnames <- colnames(data_train)[1:p]  # Extract covariate names
  
  # Fit Cox model for survival times
  fmla_T <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
  mdl_T <- coxph(fmla_T, data = data_train, x = TRUE, method = 'breslow')
  
  # Fit Cox model for censoring times
  fmla_C <- as.formula(paste("Surv(censored_T, event_C) ~ ", paste(xnames, collapse= "+")))
  mdl_C <- coxph(fmla_C, data = data_train, x = TRUE, method = 'breslow')
  
  # Determine prediction time points
  pred_time_T <- sort(unique(data_calib$censored_T[data_calib$event & 
                      data_calib$censored_T>=min(data_train$censored_T[data_train$event]) &
                      data_calib$censored_T<=max(data_train$censored_T[data_train$event])]))
  pred_time_T_test <- sort(unique(data_test$censored_T[data_test$event &
                           data_test$censored_T>=min(data_train$censored_T[data_train$event]) &
                           data_test$censored_T<=max(data_train$censored_T[data_train$event])]))
  pred_time_C <- sort(unique(data_calib$censored_T[data_calib$event_C &
                             data_calib$censored_T>=min(data_train$censored_T[data_train$event_C]) &
                             data_calib$censored_T<=max(data_train$censored_T[data_train$event_C])]))
  
  # Generate survival probability predictions
  list( 
    pred_time_T = pred_time_T,
    pred_surv_T = predictSurvProb(mdl_T, newdata = data_calib, times = pred_time_T),
    pred_time_T_test = pred_time_T_test,
    pred_surv_T_test = predictSurvProb(mdl_T, newdata = data_test, times = pred_time_T_test),
    pred_time_C = pred_time_C,
    pred_surv_C = predictSurvProb(mdl_C, newdata = data_calib, times = pred_time_C)
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
#' @return A list containing:
#'   \item{pred_time_T}{Prediction time points for T}
#'   \item{pred_surv_T}{Predicted survival probabilities for T}
#'   \item{pred_time_T_test}{Prediction time points for T (test set)}
#'   \item{pred_surv_T_test}{Predicted survival probabilities for T (test set)}
#'   \item{pred_time_C}{Prediction time points for censoring}
#'   \item{pred_surv_C}{Predicted survival probabilities for censoring}
#'
fit_models_rfsrc <- function(
    data_train, data_calib, data_test
) {
  p <- ncol(data_train) - 3  # Number of covariates
  xnames <- colnames(data_train)[1:p]  # Extract covariate names
  
  # Fit Random Survival Forest (RSF) model for survival times
  fmla_T <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
  mdl_T <- rfsrc(fmla_T, data = data_train, ntime = NULL)
  
  # Fit RSF model for censoring times
  fmla_C <- as.formula(paste("Surv(censored_T, event_C) ~ ", paste(xnames, collapse= "+")))
  mdl_C <- rfsrc(fmla_C, data = data_train, ntime = NULL)
  
  # Generate predictions for calibration and test sets
  pred_T <- predict(mdl_T, data_calib)
  pred_T_test <- predict(mdl_T, data_test)
  pred_C <- predict(mdl_C, data_calib)
  
  list( 
    pred_time_T = pred_T$time.interest,
    pred_surv_T = pred_T$survival,
    pred_time_T_test = pred_T_test$time.interest,
    pred_surv_T_test = pred_T_test$survival,
    pred_time_C = pred_C$time.interest,
    pred_surv_C = pred_C$survival
  )
}

#' Fit Super Learner survival models
#'
#' Uses Super Learner (SL) to fit an ensemble model for survival prediction.
#'
#' @param data_train Training dataset
#' @param data_calib Calibration dataset
#' @param data_test Test dataset
#' @return A list containing:
#'   \item{pred_time_T}{Prediction time points for T}
#'   \item{pred_surv_T}{Predicted survival probabilities for T}
#'   \item{pred_time_T_test}{Prediction time points for T (test set)}
#'   \item{pred_surv_T_test}{Predicted survival probabilities for T (test set)}
#'   \item{pred_time_C}{Prediction time points for censoring}
#'   \item{pred_surv_C}{Predicted survival probabilities for censoring}
#'
fit_models_sl <- function(
    data_train, data_calib, data_test
) {
  p <- ncol(data_train) - 3  # Number of covariates
  
  # Define the library of models for Super Learner
  SL.library <- c("survSL.km", "survSL.coxph", "survSL.expreg", 
                  "survSL.weibreg", "survSL.loglogreg", "survSL.gam", 
                  "survSL.rfsrc")
  
  # Fit Super Learner model for survival
  fit <- survSuperLearner(
    time = data_train$censored_T,
    event = data_train$event,
    X = data_train[, 1:p],
    newX = data_calib[, 1:p],
    new.times = sort(unique(data_calib$censored_T[data_calib$event])),
    event.SL.library = SL.library,          
    cens.SL.library = SL.library,
    verbose = FALSE
  )
  
  # Generate predictions for test and calibration sets
  pred_time_T_test <- sort(unique(data_test$censored_T[data_test$event]))
  pred_T_test <- predict.survSuperLearner(fit, newdata = data_test[, 1:p],
                                          new.times = pred_time_T_test,
                                          onlySL = TRUE)
  pred_time_C <- sort(unique(data_calib$censored_T[!data_calib$event]))
  pred_C <- predict.survSuperLearner(fit, newdata = data_calib[, 1:p],
                                     new.times = pred_time_C,
                                     onlySL = TRUE)
  
  list( 
    pred_time_T = sort(unique(data_calib$censored_T[data_calib$event])),
    pred_surv_T = fit$event.SL.predict,
    pred_time_T_test = pred_time_T_test,
    pred_surv_T_test = pred_T_test$event.SL.predict,
    pred_time_C = pred_time_C,
    pred_surv_C = pred_C$cens.SL.predict
  )
}

