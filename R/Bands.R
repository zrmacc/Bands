# Purpose: Master banding function
# Updated: 19/03/20

#' Confidence Band for Survival Functions
#' 
#' @param time Observation time.
#' @param status Status, 1 for an event, 0 for censoring.
#' @param arm Treatment indicator, for two-sample problem: 1 for target group, 0
#'   for reference.
#' @param sig Significance level.
#' @param paths Sample paths. 
#' @param weights Either equiprecision "ep" or Hall-Wellner "hw". 
#' @param L Lower endpoint of interval over which confidence bands are sought. 
#' @param U Upper endpoint of interval over which confidence bands are sought. 
#' 
#' @return Object of class band, containing a table with the lower and upper
#'   confidence bands.
#' @export

survBands = function(time,status=NULL,arm=NULL,sig=0.05,paths=1e3,weights="ep",L=NULL,U=NULL){
  ## Input checks
  # Positivity
  if(min(time)<0){
    stop("Non-negative observation times required.");
  };
  # Missingness
  if(sum(is.na(time))+sum(is.na(status))>0){
    stop("Input contains missing data");
  };
  # Weights
  if(!(weights %in% c("ep","hw"))){
    stop("Select weights from 'ep' or 'hw'");
  };
  
  # Status
  n = length(time);
  if(is.null(status)){status=rep(1,n)};
  status.levels = sort(unique(status));
  if(length(status.levels)==1){
    if(status.levels!=1){status=rep(1,n)};
  }
  if(length(status.levels)==2){
    if(!all.equal(status.levels,c(0,1))){stop("Numeric 0,1 coding is expected for status.")};
  }
  if(length(status.levels)>2){stop("Only two levels are expected for status.")};
  
  # Arm
  if(!is.null(arm)){
    arm.levels = sort(unique(arm));
    if(length(arm.levels)==1){
      arm = NULL;
    };
    if(length(arm.levels)==2){
      if(!all.equal(arm.levels,c(0,1))){stop("Numeric 0,1 coding is expected for arm.")};
    }
  };
  
  ## Fit confidence bands
  Out = NULL;
  # Arm
  if(is.null(arm)){
    Out = oneSampleBand(time=time,status=status,sig=sig,paths=paths,weights=weights,L=L,U=U);
  } else {
    Out = twoSampleBand(time=time,status=status,arm=arm,sig=sig,paths=paths,weights=weights,L=L,U=U)
  }
  return(Out);
}
