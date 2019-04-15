#' Fitted Survival Bands
#'
#' Defines the object class returned by fitting functions.
#'
#' @slot Alpha Alpha level. 
#' @slot Crit Critical value.
#' @slot Paths Sample paths. 
#' @slot Samples Number of samples. 
#' @slot Table Results table. 
#' @name band-class
#' @rdname band-class
#' @exportClass band

setClass(Class="band",representation=representation(Alpha="numeric",Crit="numeric",Paths="numeric",Samples="numeric",Table="data.frame"));

########################
# Plot Method
########################

#' Plot Method for Survival Band
#'
#' Plot method for objects of class \code{band}.
#'
#' @param x An object of class \code{band}.
#' @param y Unused.
#' @param ... Unused.
#' @export

plot.band = function(x,y,...){
  q = plotBands(x);
  return(q);
}

########################
# Print Method
########################

#' Print Method for Survival Band
#'
#' Print method for objects of class \code{band}.
#'
#' @param x An object of class \code{band}.
#' @param ... Unused.
#' @export

print.band = function(x,...){
  # Samples
  cat(x@Samples,"Sample Survival Band.\n");
  # Alpha level
  cat("Significance level:",x@Alpha,"\n");
  # Critical value
  cat("Critical value:",round(x@Crit,digits=2),"\n");
  # Sample paths
  cat("Sample paths:",x@Paths,"\n");
}

########################
# Show Method
########################

#' Show Method for Survival Band
#'
#' @param object An object of class \code{band}.
#' @rdname band-method
#' @importFrom methods show

setMethod(f="show",signature=c(object="band"),definition=function(object){print.band(x=object)});