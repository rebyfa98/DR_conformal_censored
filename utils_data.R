########################################
## utils_data.R
##
## Generates input survival data in the format needed for the method.
########################################

#' Simulate exponential survival data
#'
#' @param seed Random seed
#' @param n Sample size
#' @param beta_T Coefficients for generating survival time T
#' @param beta_C Coefficients for generating survival time C
#' @param lambda_T Rate for generating survival time T
#' @param lambda_C Rate for generating censoring time C
#' @return A list containing:
#'   \item{data}{Dataset containing covariates, censored survival times, 
#'   and event indicators(data frame)}
#'   \item{T_true}{True survival times before censoring (numeric vector)}
#'
simulate_exponential_data <- function(
    seed,
    n = 3000,         # Total sample size
    beta_T = c(-1, 1),# Regression coefficients for survival time generation
    beta_C = c(0, 0), # Regression coefficients for censoring time generation
    lambda_T = 0.01,  # Rate parameter for survival time distribution
    lambda_C = 0.01  # Rate parameter for censoring time distribution
) {
  set.seed(seed)
  p <- length(beta_T) # Number of covariates
  
  # Helpers for data simulation
  gen_X <- function(p, n) rnorm(p * n, 0, 1)
  gen_T <- function(lambda, x, beta) rexp(1, lambda * exp(x %*% beta))
  gen_C <- function(lambda, x, beta) rexp(1, lambda * exp(x %*% beta))
  
  # Generate data
  X <- matrix(gen_X(p, n), ncol = p)
  T <- apply(X, 1, gen_T, lambda = lambda_T, beta = beta_T)
  C <- apply(X, 1, gen_C, lambda = lambda_C, beta = beta_C)
  
  # Determine censoring status
  event    <- (T <= C)
  event_C  <- (T > C)
  censored_T <- pmin(T, C)
  
  # Organize into a data.frame
  data <- data.frame(X = X, 
                     censored_T = censored_T, 
                     event = event, 
                     event_C = event_C)
  
  # Name covariates
  xnames <- paste0("X", 1:p)
  colnames(data)[1:p] <- xnames
  
  list(
    data = data,
    T_true = T # Keep track of true survival times
  )
}

#' Simulate lognormal survival data
#'
#' Generates synthetic survival data where survival and censoring times follow a lognormal distribution.
#'
#' @param seed Random seed for reproducibility
#' @param n Sample size
#' @param p Number of covariates
#' @param A_T A function defining the treatment effect condition based on covariates `x`.
#' @param A_C A function defining the censoring condition based on covariates `x`.
#' @return A list containing:
#'   \item{data}{A data frame with generated covariates, censored survival times, and event indicators.}
#'   \item{T_true}{A numeric vector of true survival times before censoring.}
#'
simulate_lognorm_data <- function(
    seed,
    n = 3000,
    p = 100,
    A_T = function(x) (x[2] < 0 & x[3] > 0 & x[4] > 0),
    A_C = function(x) (x[1] > 0)
) {
  set.seed(seed)
  
  # Helpers for data simulation
  gen_X <- function(p, n) runif(p * n, -1, 1)  # Generate independent covariates
  
  gen_T <- function(x) {
    if (A_T(x)) {
      return(exp(rnorm(1, log(10), 1)))  # Shorter survival time
    } else {
      return(exp(rnorm(1, log(1000), 1)))  # Longer survival time
    }
  }
  
  gen_C <- function(x) {
    if (A_C(x)) {
      return(exp(rnorm(1, log(10), 1)))  # Shorter censoring time
    } else {
      return(exp(rnorm(1, log(1000), 1)))  # Longer censoring time
    }
  }
  
  # Generate data
  X <- matrix(gen_X(p, n), ncol = p)
  T <- apply(X, 1, gen_T)  # True survival times
  C <- apply(X, 1, gen_C)  # Censoring times
  
  # Determine censoring status
  event    <- (T <= C)  # Event indicator (1 = observed, 0 = censored)
  event_C  <- (T > C)   # Censoring indicator (opposite of event)
  censored_T <- pmin(T, C)  # Observed survival times (either T or C)
  
  # Organize into a data frame
  data <- data.frame(X = X, 
                     censored_T = censored_T, 
                     event = event, 
                     event_C = event_C)
  
  # Name covariates
  xnames <- paste0("X", 1:p)
  colnames(data)[1:p] <- xnames
  
  list(
    data = data,
    T_true = T  # Keep track of true survival times before censoring
  )
}
