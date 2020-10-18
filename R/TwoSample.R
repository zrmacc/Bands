# Purpose: Estimates 2-sample confidence bands
# Updated: 2020-10-17

#' Estimate 2-Sample Confidence Bands
#'
#' @param time Observation time.
#' @param status Status, 1 for an event, 0 for censoring.
#' @param arm Treatment indicator, 1 for target group, 0 for reference.
#' @param alpha Type I error level.
#' @param paths Sample paths.
#' @param method Either equiprecision "EP" or Hall-Wellner "HW".
#' @param lower Lower endpoint of interval over which confidence bands are sought.
#' @param upper Upper endpoint of interval over which confidence bands are sought.
#'
#' @return Object of class band, containing the lower and upper confidence bands.
#'
#' @importFrom data.table data.table
#' @importFrom methods new
#' @importFrom stats quantile rnorm stepfun
#' @importFrom survival Surv survfit

Bands.TwoSample <- function(
  time, 
  status, 
  arm, 
  alpha = 0.05, 
  paths = 1e3, 
  method = "EP", 
  lower = NULL, 
  upper = NULL
) {
  
  # Sample size
  n <- length(time)
  
  # KM Curves
  km_0 <- survfit(Surv(time[arm == 0], status[arm == 0]) ~ 1)
  km_1 <- survfit(Surv(time[arm == 1], status[arm == 1]) ~ 1)
  
  # Lower endpoint, no less than the max of the minimum observed event times.
  lower_min <- max(min(km_0$time[km_0$surv < 1]), min(km_1$time[km_1$surv < 1]))
  if (is.null(lower) || lower < lower_min) {
    lower <- lower_min
  }
  
  # Upper endpoint
  upper_max <- min(max(km_0$time[km_0$surv > 0]), max(km_1$time[km_1$surv > 0]))
  if (is.null(upper) || upper > upper_max) {
    upper <- upper_max
  }
  
  # Jump indicators
  is_jump_0 <- (km_0$n.event > 0) & (km_0$time >= lower) & (km_0$time <= upper)
  is_jump_1 <- (km_1$n.event > 0) & (km_1$time >= lower) & (km_1$time <= upper)
  
  # Event times.
  tab_0 <- data.table("time" = unique(km_0$time[is_jump_0]), key = "time")
  tab_1 <- data.table("time" = unique(km_1$time[is_jump_1]), key = "time")
  
  # Events counts.
  tab_0$events0 <- km_0$n.event[is_jump_0]
  tab_1$events1 <- km_1$n.event[is_jump_1]
  process_tab <- merge(x = tab_0, y = tab_1, all = TRUE)
  process_tab[is.na(process_tab)] <- 0
  
  # NARs.
  nar_0 <- stepfun(x = km_0$time, y = c(km_0$n, km_0$n.risk))
  nar_1 <- stepfun(x = km_1$time, y = c(km_1$n, km_1$n.risk))
  process_tab$nar0 <- nar_0(process_tab$time)
  process_tab$nar1 <- nar_1(process_tab$time)
  
  # Survival rates.
  surv_0 <- stepfun(x = km_0$time, y = c(1, km_0$surv))
  surv_1 <- stepfun(x = km_1$time, y = c(1, km_1$surv))
  process_tab$surv0 <- surv_0(process_tab$time)
  process_tab$surv1 <- surv_1(process_tab$time)
  
  # Add cumulative hazard to table. 
  process_tab$cumhaz0 <- cumsum(process_tab$events0 / process_tab$nar0)
  process_tab$cumhaz1 <- cumsum(process_tab$events1 / process_tab$nar1)
  
  # SE step functions.
  se_0 <- stepfun(x = km_0$time, y = c(0, km_0$std.err))
  se_1 <- stepfun(x = km_1$time, y = c(1, km_1$std.err))
  process_tab$se0 <- se_0(process_tab$time)
  process_tab$se1 <- se_1(process_tab$time)
  
  # Weights.
  if (method == "EP") {
    process_tab$weight <- 1
  } else if (method == "HW") {
    process_tab$weight <- 1 / (1 + process_tab$se0^2 + process_tab$se1^2)
  }
  
  # Observed KS statistic. 
  process_tab$delta <- process_tab$surv1 - process_tab$surv0
  ks_obs <- sqrt(n) * max(process_tab$weight * abs(process_tab$delta))
  
  # Simulations:
  # Generate sample paths and find the supremum.
  sim <- sapply(
    seq_len(paths),
    function(x) {
      
      # Arm 0.
      sample_path_0 <- SamplePath(
        times = process_tab$time0,
        events = process_tab$events0,
        nar = process_tab$nar0,
        surv = process_tab$surv0,
        weight = process_tab$weight
      )
      
      # Arm 1.
      sample_path_1 <- SamplePath(
        times = process_tab$time1,
        events = process_tab$events1,
        nar = process_tab$nar1,
        surv = process_tab$surv1,
        weight = process_tab$weight
      )
      
      # Difference.
      sample_path <- sample_path_1 - sample_path_0
      sup <- sqrt(n) * max(abs(sample_path))
      return(sup)
    }
  )

  # Critical value.
  critical_value <- as.numeric(quantile(sim, probs = 1 - alpha))
  
  # Confidence bands.
  process_tab$lower <- process_tab$delta - critical_value / (process_tab$weight * sqrt(n))
  process_tab$upper <- process_tab$delta + critical_value / (process_tab$weight * sqrt(n))
  
  # P-value.
  p_val <- mean(c(1, sim >= ks_obs))
  
  # Output
  out <- new(
    Class = "band", 
    Alpha = alpha, 
    Crit = critical_value, 
    Paths = paths, 
    Pvalue = p_val,
    Samples = 2, 
    Table = process_tab
  )
  return(out)
}
