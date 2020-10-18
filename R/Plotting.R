# Purpose: Plotting Confidence Bands
# Updated: 19/03/20

#' Master Plotting Function
#'
#' @param band Object of class band.
#' @return An object of class `ggplot`.

PlotBands <- function(band) {
  samples <- band@Samples
  q <- NULL
  if (samples == 1) {
    q <- Plot.OneSample(band)
  } else if (samples == 2) {
    q <- Plot.TwoSample(band)
  }
  return(q)
}

#' Plot 1-Sample Confidence Band
#'
#' @param band Object of class band.
#' @import ggplot2
#' @return An object of class `ggplot`.

Plot.OneSample <- function(band) {
  q <- ggplot(data = band@Table) +
    theme_bw() + 
    geom_ribbon(
      aes(x = time, ymin = lower, ymax = upper), 
      fill = "#0073C2FF", 
      alpha = 0.5
    ) + geom_step(
      aes(x = time, y = surv), 
      color = "#0073C2FF"
    ) +
    xlab("Time") + 
    ylab("Survival")
  return(q)
}

#' Plot 2-Sample Confidence Band
#'
#' @param B Object of class band.
#' @import ggplot2
#' @return An object of class `ggplot`.

Plot.TwoSample <- function(band) {
  q <- ggplot(data = band@Table) +
    theme_bw() + 
    geom_step(
      aes(x = time, y = delta), 
      color = "#0073C2FF"
    ) + geom_ribbon(
      aes(x = time, ymin = lower, ymax = upper), 
      fill = "#0073C2FF", 
      alpha = 0.5
    ) + geom_hline(
      yintercept = 0, 
      linetype = "dashed", 
      color = "gray"
    ) + 
    xlab("Time") + 
    ylab(expression(Survival ~ Difference ~ (S[1] - S[0])))
  return(q)
}
