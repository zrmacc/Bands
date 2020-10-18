#' Fitted Survival Bands
#'
#' Defines the object class returned by fitting functions.
#'
#' @slot Alpha Alpha level.
#' @slot Crit Critical value.
#' @slot Paths Number of simulated sample paths.
#' @slot Pvalue P-value for the difference of survival curves.
#' @slot Samples Number of samples.
#' @slot Table Tabulated survival function with confidence bands.
#' @name band-class
#' @rdname band-class
#' @exportClass band

setClass(
  Class = "band", 
  representation = representation(
    Alpha = "numeric", 
    Crit = "numeric", 
    Paths = "numeric", 
    Pvalue = "numeric",
    Samples = "numeric", 
    Table = "data.frame"
  )
)


# -----------------------------------------------------------------------------
# Plot Method
# -----------------------------------------------------------------------------

#' Plot Method for Survival Band
#'
#' Plot method for objects of class \code{band}.
#'
#' @param x An object of class \code{band}.
#' @param y Unused.
#' @param ... Unused.
#' @export

plot.band <- function(x, y, ...) {
  q <- PlotBands(x)
  return(q)
}

# -----------------------------------------------------------------------------
# Print Method
# -----------------------------------------------------------------------------

#' Print Method for Survival Band
#'
#' Print method for objects of class \code{band}.
#'
#' @param x An object of class \code{band}.
#' @param ... Unused.
#' @export

print.band <- function(x, ...) {
  
  # Samples
  samples <- x@Samples
  cat(samples, "Sample Survival Band.\n")
  
  # Alpha level
  cat("Significance level:", x@Alpha, "\n")
  
  # Critical value
  cat("Critical value:", round(x@Crit, digits = 3), "\n")
  
  # Sample paths
  cat("Sample paths:", x@Paths, "\n")
  
  # P-value
  cat("P-value for difference of survival curves:", 
      signif(x@Pvalue, digits = 3), "\n")
}


# -----------------------------------------------------------------------------
# Show Method
# -----------------------------------------------------------------------------

#' Show Method for Survival Band
#'
#' @param object An object of class \code{band}.
#' @rdname band-method
#' @importFrom methods show

setMethod(
  f = "show", 
  signature = c(object = "band"), 
  definition = function(object) {
    print.band(x = object)
})
