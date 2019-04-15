# Purpose: Estimates 1-sample confidence bands

#' Estimate 1-Sample Confidence Bands
#' 
#' Constructs the confidence band on the log cumulative hazard (CH) scale,
#' then transforms onto the survival scale. The critical value, on the
#' log CH scale, is determined via simulation. 
#' 
#' @param time Observation time.
#' @param status Status, 1 for an event, 0 for censoring.
#' @param sig Significance level.
#' @param paths Sample paths. 
#' @param weights Either equiprecision "ep" or Hall-Wellner "hw". 
#' @param L Lower endpoint of interval over which confidence bands are sought. 
#' @param U Upper endpoint of interval over which confidence bands are sought.  
#' 
#' @return Object of class band, containing a table with the lower and upper
#'   confidence bands.
#' 
#' @importFrom data.table data.table
#' @importFrom methods new
#' @importFrom stats quantile rnorm
#' @importFrom survival Surv survfit 

oneSampleBand = function(time,status,sig=0.05,paths=1e3,weights="ep",L=NULL,U=NULL){
  # Sample size
  n = length(time);
  # Observed KM
  K = survfit(Surv(time,status)~1);
  
  # Lower endpoint
  L.min = min(K$time[K$surv<1]);
  if(is.null(L)){L=L.min};
  if(L<L.min){L=L.min};
  # Upper endpoint
  U.max = max(K$time[K$surv>0]);
  if(is.null(U)){U=U.max};
  if(U>U.max){U=U.max};
  
  # Key
  flag = (K$n.event>0)&(K$time>=L)&(K$time<=U);
  m = sum(flag);
  
  # Times
  Tab = data.table("time"=unique(K$time[flag]));
  
  # Events
  Tab$events = K$n.event[flag];

  # NARs
  Tab$nar = K$n.risk[flag];
  
  # Survival estimates
  Tab$surv = K$surv[flag];
  
  # Cumulative hazard
  Tab$cumhaz = cumsum(Tab$events/Tab$nar);
  
  # Standard errors
  Tab$se = K$std.err[flag];
  
  # Weights
  if(weights=="ep"){
    Tab$weight = 1;
  } else if(weights=="hw"){
    Tab$weight = 1/(1+Tab$se^2);
  }

  # "Loop"
  aux = function(r){
    # Weights
    if(max(Tab$events)==1){
      Z = rnorm(n=m);
    } else {
      Z = sapply(E,function(x){rnorm(n=1,mean=0,sd=sqrt(x))});
    }
    # Loop over doots
    V = (Z/Tab$nar);
    # Cumulative sum
    Q = cumsum(V);
    # Scale
    R = -sqrt(n)*(Tab$surv)*(Tab$weight)*Q;
    # Output
    return(max(abs(R)))
  }
  
  # Simulations
  Sups = sapply(seq(1:paths),aux);
  Ka = as.numeric(quantile(Sups,probs=1-sig));
  
  # Results
  Tab$L = sapply(Tab$surv-Ka/(Tab$weight*sqrt(n)),function(x){max(x,0)});
  Tab$U = sapply(Tab$surv+Ka/(Tab$weight*sqrt(n)),function(x){min(x,1)});
  
  # Output
  Out = new(Class="band",Alpha=sig,Crit=Ka,Paths=paths,Samples=1,Table=Tab);
  return(Out);
}
