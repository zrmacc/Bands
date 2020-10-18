# Purpose: Estimates 1-sample confidence bands
# Updated: 2020-10-17

#' Sample Path of Survival Process
#' 
#' Perturbation realization of a sample path from \eqn{\sqrt{n}w(t)\{\hat{S}(t)-S(t)\}}.
#' 
#' @param times Distinct observed event times.
#' @param events Number of events at each time.
#' @param nar Number at risk.
#' @param surv Survival probability. 
#' @param weight Weight to give each time point.
#' @importFrom stats rnorm
#' @return Numeric vector of the same length as `times`.

SamplePath <- function(
  times,
  events,
  nar,
  surv,
  weight = NULL
) {
  
  # Default to uniform weights.
  if (is.null(weight)) {
    weight <- rep(1, length(times))
  }
  
  # Perturbation weights.
  z <- sapply(events, function(x) {
    rnorm(n = 1, mean = 0, sd = sqrt(x))
  })
  
  # Scale
  sample_path <- -weight * surv * cumsum(z * events / nar)
  
  # Output
  return(sample_path)
}

#' Estimate 1-Sample Confidence Bands
#'
#' Constructs the confidence band on the log cumulative hazard (CH) scale,
#' then transforms onto the survival scale. The critical value, on the
#' log CH scale, is determined via simulation.
#'
#' @param time Observation time.
#' @param status Status, 1 for an event, 0 for censoring.
#' @param alpha Alpha level.
#' @param paths Sample paths.
#' @param method Either equiprecision "EP" or Hall-Wellner "HW".
#' @param lower Lower endpoint of interval over which confidence bands are sought.
#' @param upper Upper endpoint of interval over which confidence bands are sought.
#'
#' @return Object of class band, containing a table with the lower and upper
#'   confidence bands.
#'
#' @importFrom data.table data.table
#' @importFrom methods new
#' @importFrom stats quantile
#' @importFrom survival Surv survfit

Bands.OneSample <- function(
  time, 
  status, 
  alpha = 0.05, 
  paths = 1e3, 
  method = "EP", 
  lower = NULL, 
  upper = NULL
) {
  
  # Sample size
  n <- length(time)
  
  # Observed KM
  km_fit <- survfit(Surv(time, status) ~ 1)

  # Lower endpoint, no less than the minimum observed event time.
  lower_min <- min(km_fit$time[km_fit$surv < 1])
  if (is.null(lower) || lower < lower_min) {
    lower <- lower_min
  }
  
  # Upper endpoint, no greater than the maximum observed event time.
  upper_max <- max(km_fit$time[km_fit$surv > 0])
  if (is.null(upper) || upper > upper_max) {
    upper <- upper_max
  }

  # Jump indicator.
  is_jump <- (km_fit$n.event > 0) & (km_fit$time >= lower) & (km_fit$time <= upper)
  
  # Process table.
  process_tab <- data.table(
    "time" = unique(km_fit$time[is_jump]),
    "events" = km_fit$n.event[is_jump],
    "nar" = km_fit$n.risk[is_jump],
    "surv" = km_fit$surv[is_jump],
    "cumhaz" = km_fit$cumhaz[is_jump],
    "se" = km_fit$std.err[is_jump]
  )

  # Weights
  if (method == "EP") {
    process_tab$weight <- 1
  } else if (method == "HW") {
    process_tab$weight <- 1 / (1 + process_tab$se^2)
  }

  # Simulations:
  # Generate a sample path and find the supremum.
  sim <- sapply(
    seq_len(paths),
    function(x) {
      sample_path <- SamplePath(
        times = process_tab$time,
        events = process_tab$events,
        nar = process_tab$nar,
        surv = process_tab$surv,
        weight = process_tab$weight
      )
      sup <- sqrt(n) * max(abs(sample_path))
      return(sup)
    }
  )
  
  # Critical value.
  critical_value <- as.numeric(quantile(sim, probs = 1 - alpha))

  # Confidence bands.
  process_tab$lower <- sapply(
    process_tab$surv - critical_value / (process_tab$weight * sqrt(n)),
    function(x) {max(x, 0)}
  )
  process_tab$upper <- sapply(
    process_tab$surv + critical_value / (process_tab$weight * sqrt(n)),
    function(x) {min(x, 1)}
  )

  # Output
  out <- new(
    Class = "band", 
    Alpha = alpha, 
    Crit = critical_value, 
    Paths = paths, 
    Pvalue = NA_real_,
    Samples = 1, 
    Table = process_tab
  )
  return(out)
}
