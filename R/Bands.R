# Purpose: Overall banding function.
# Updated: 2020-10-17

#' Confidence Band for Survival Curves
#'
#' @param time Observation time.
#' @param status Status, 1 for an event, 0 for censoring.
#' @param arm Treatment indicator, for two-sample problem: 1 for target group, 0
#'   for reference.
#' @param alpha Type 1 error level.
#' @param paths Sample paths.
#' @param method Either equiprecision "EP" or Hall-Wellner "HW".
#' @param lower Lower endpoint of interval over which confidence bands are sought.
#' @param upper Upper endpoint of interval over which confidence bands are sought.
#'
#' @return Object of class band, containing a table with the lower and upper
#'   confidence bands.
#' @export

SurvBands <- function(
  time, 
  status = NULL, 
  arm = NULL, 
  alpha = 0.05, 
  paths = 1e3, 
  method = "EP", 
  lower = NULL, 
  upper = NULL) {
  
  ## Input checks
  # Missingness
  if (sum(is.na(time)) + sum(is.na(status)) > 0) {
    stop("Input contains missing data")
  }
  
  # Method
  if (!(method %in% c("EP", "HW"))) {
    stop("Select method from 'EP' or 'HW'")
  }

  # Ensure status is coded 0/1.
  n <- length(time)
  if (is.null(status)) {
    status <- rep(1, n)
  }
  status_levels <- sort(unique(status))
  if (length(status_levels) == 1) {
    if (status_levels != 1) {
      status <- rep(1, n)
    }
  }
  if (length(status_levels) == 2) {
    if (!all.equal(status_levels, c(0, 1))) {
      stop("Numeric 0,1 coding is expected for status.")
    }
  }
  if (length(status_levels) > 2) {
    stop("Only two levels are expected for status.")
  }

  # Ensure arm is coded 0/1.
  if (!is.null(arm)) {
    arm_levels <- sort(unique(arm))
    if (length(arm_levels) == 1) {
      arm <- NULL
    }
    if (length(arm_levels) == 2) {
      if (!all.equal(arm_levels, c(0, 1))) {
        stop("Numeric 0,1 coding is expected for arm.")
      }
    }
  }

  # Confidence bands.
  out <- NULL
  if (is.null(arm)) {
    out <- Bands.OneSample(
      time = time,
      status = status,
      alpha = alpha,
      paths = paths,
      method = method,
      lower = lower, 
      upper = upper
    ) 
  } else {
    out <- Bands.TwoSample(
      time = time,
      status = status,
      arm = arm,
      alpha = alpha,
      paths = paths,
      method = method,
      lower = lower, 
      upper = upper
    ) 
  }
  return(out)
}
