########################################
## utils_data.R
##
## Generates input survival data in the format needed for the method.
########################################

#' Build a survival dataset in the common format used by the method
#'
#' @param X Covariate matrix
#' @param T_true True event times
#' @param C_true True censoring times
#' @return A list containing:
#'   \item{data}{Dataset containing covariates, censored survival times,
#'   and event indicators (data frame)}
#'   \item{T_true}{True survival times before censoring (numeric vector)}
#'
build_survival_dataset <- function(X, T_true, C_true) {
  # Determine censoring status
  event <- T_true <= C_true

  # Organize into a data.frame
  data <- data.frame(
    X = X,
    censored_T = pmin(T_true, C_true),
    event = event,
    event_C = !event
  )

  # Name covariates
  p <- ncol(X)
  colnames(data)[seq_len(p)] <- paste0("X", seq_len(p))

  list(
    data = data,
    T_true = T_true
  )
}

#' Simulate survival data for setting 1
#'
#' @param seed Random seed
#' @param n Sample size
#' @return A list containing the simulated dataset and the true event times
#'
simulate_setting1 <- function(seed, n = 3000) {
  set.seed(seed)

  p <- 2

  # Generate data
  X <- matrix(rnorm(n * p), ncol = p)
  T_true <- apply(X, 1, function(x) rexp(1, rate = exp(-x[1] + x[2])))
  C_true <- rexp(n, rate = 1 / 3)

  build_survival_dataset(X, T_true, C_true)
}

#' Simulate survival data for setting 2
#'
#' @param seed Random seed
#' @param n Sample size
#' @return A list containing the simulated dataset and the true event times
#'
simulate_setting2 <- function(seed, n = 3000) {
  set.seed(seed)

  p <- 10

  # Generate data
  X <- matrix(rnorm(n * p), ncol = p)
  T_true <- apply(
    X,
    1,
    function(x) rexp(1, rate = (2 / 3) * exp((x[1] * x[2] - x[3]^2) / 3))
  )
  C_true <- apply(
    X,
    1,
    function(x) rexp(1, rate = (2 / 3) * exp((x[3] - x[4]^2) / 3))
  )

  build_survival_dataset(X, T_true, C_true)
}

#' Simulate survival data for setting 3
#'
#' @param seed Random seed
#' @param n Sample size
#' @return A list containing the simulated dataset and the true event times
#'
simulate_setting3 <- function(seed, n = 3000) {
  set.seed(seed)

  p <- 100

  # Generate data
  X <- matrix(runif(n * p, min = -1, max = 1), ncol = p)

  A_T <- function(x) all(x[1:5] > 0) && all(x[6:10] < 0)
  A_C <- function(x) x[1] > 0 && x[2] < 0

  T_true <- apply(
    X,
    1,
    function(x) exp(rnorm(1, mean = if (A_T(x)) log(10) else log(1000), sd = 1))
  )
  C_true <- apply(
    X,
    1,
    function(x) exp(rnorm(1, mean = if (A_C(x)) log(10) else log(1000), sd = 1))
  )

  build_survival_dataset(X, T_true, C_true)
}

#' Simulate survival data for setting 4
#'
#' @param seed Random seed
#' @param n Sample size
#' @return A list containing the simulated dataset and the true event times
#'
simulate_setting4 <- function(seed, n = 3000) {
  set.seed(seed)

  p <- 100

  # Generate data
  X <- matrix(runif(n * p, min = -1, max = 1), ncol = p)

  A <- function(x) x[2] < 0 && x[3] > 0 && x[4] > 0
  A_C <- function(x) x[1] < 0

  T_true <- apply(
    X,
    1,
    function(x) exp(rnorm(1, mean = if (A(x)) log(10) else log(1000), sd = 1))
  )
  C_true <- apply(
    X,
    1,
    function(x) exp(rnorm(1, mean = if (A_C(x)) log(10) else log(1000), sd = 1))
  )

  build_survival_dataset(X, T_true, C_true)
}

#' Simulate survival data for setting 5
#'
#' @param seed Random seed
#' @param n Sample size
#' @return A list containing the simulated dataset and the true event times
#'
simulate_setting5 <- function(seed, n = 3000) {
  set.seed(seed)

  p <- 100

  # Generate data
  X <- matrix(runif(n * p), ncol = p)

  T_true <- apply(
    X,
    1,
    function(x) {
      meanlog <- (x[1] - 0.5)^2 +
        x[2] * x[3] -
        as.numeric(x[3] < 0.5 && x[4] > 0.5) +
        sqrt(x[5]) +
        (x[6] + x[7] - 0.5)^3
      rlnorm(1, meanlog = meanlog, sdlog = 1)
    }
  )
  C_true <- apply(
    X,
    1,
    function(x) {
      meanlog <- (x[1] + x[2] - 1)^2 -
        x[3] * x[4] +
        as.numeric(x[6] > 0.5) -
        (x[7] - 0.5)^3 * x[8]
      rlnorm(1, meanlog = meanlog, sdlog = 1)
    }
  )

  build_survival_dataset(X, T_true, C_true)
}

#' Simulate survival data for setting 6
#'
#' @param seed Random seed
#' @param n Sample size
#' @return A list containing the simulated dataset and the true event times
#'
simulate_setting6 <- function(seed, n = 3000) {
  set.seed(seed)

  p <- 100

  # Generate data
  X <- matrix(runif(n * p), ncol = p)

  T_true <- apply(
    X,
    1,
    function(x) {
      meanlog <- 0.126 * (x[1] + sqrt(x[3] * x[5])) + 1
      sdlog <- (x[2] + 2) / 4
      rlnorm(1, meanlog = meanlog, sdlog = sdlog)
    }
  )
  C_true <- apply(X, 1, function(x) rexp(1, rate = x[6] / 2))

  build_survival_dataset(X, T_true, C_true)
}

#' Simulate survival data for one of the six settings
#'
#' @param setting A string specifying the synthetic setting
#' @param seed Random seed
#' @param n Sample size
#' @return A list containing the simulated dataset and the true event times
#'
simulate_data <- function(setting, seed, n = 3000) {
  setting <- match.arg(
    setting,
    choices = paste0("setting", 1:6)
  )

  switch(
    setting,
    setting1 = simulate_setting1(seed = seed, n = n),
    setting2 = simulate_setting2(seed = seed, n = n),
    setting3 = simulate_setting3(seed = seed, n = n),
    setting4 = simulate_setting4(seed = seed, n = n),
    setting5 = simulate_setting5(seed = seed, n = n),
    setting6 = simulate_setting6(seed = seed, n = n)
  )
}
