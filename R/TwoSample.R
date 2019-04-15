# Purpose: Estimates 2-sample confidence bands
# Updated: 19/03/20

#' Estimate 2-Sample Confidence Bands
#' 
#' @param time Observation time.
#' @param status Status, 1 for an event, 0 for censoring.
#' @param arm Treatment indicator, 1 for target group, 0 for reference.
#' @param sig Significance level.
#' @param paths Sample paths. 
#' @param weights Either equiprecision "ep" or Hall-Wellner "hw". 
#' @param L Lower endpoint of interval over which confidence bands are sought. 
#' @param U Upper endpoint of interval over which confidence bands are sought. 
#' 
#' @return Object of class band, containing the lower and upper confidence bands.
#' 
#' @importFrom data.table data.table 
#' @importFrom methods new
#' @importFrom stats quantile rnorm stepfun
#' @importFrom survival Surv survfit 

twoSampleBand = function(time,status,arm,sig=0.05,paths=1e3,weights="ep",L=NULL,U=NULL){
  # Sample size
  n = length(time);
  # KM Curves
  K0 = survfit(Surv(time[arm==0],status[arm==0])~1);
  K1 = survfit(Surv(time[arm==1],status[arm==1])~1);
  
  # Lower endpoint
  L.min = max(min(K0$time[K0$surv<1]),min(K1$time[K1$surv<1]));
  if(is.null(L)){L=L.min};
  if(L<L.min){L=L.min};
  # Upper endpoint
  U.max = min(max(K0$time[K0$surv>0]),max(K1$time[K1$surv>0]));
  if(is.null(U)){U=U.max};
  if(U>U.max){U=U.max};
  
  # Keys
  flag0 = (K0$n.event>0)&(K0$time>=L)&(K0$time<=U);
  flag1 = (K1$n.event>0)&(K1$time>=L)&(K1$time<=U);
  
  # Event times
  E0 = data.table("time"=unique(K0$time[flag0]),key="time");
  E1 = data.table("time"=unique(K1$time[flag1]),key="time");
  
  # Events
  E0$events0 = K0$n.event[flag0];
  E1$events1 = K1$n.event[flag1];
  
  # Merge
  Tab = merge(x=E0,y=E1,all=T);
  Tab[is.na(Tab)] = 0;
  rm(E0,E1);

  # NAR step functions
  f0 = stepfun(x=K0$time,y=c(K0$n,K0$n.risk));
  f1 = stepfun(x=K1$time,y=c(K1$n,K1$n.risk));
  
  # Add NARs to table
  Tab$nar0 = f0(Tab$time);
  Tab$nar1 = f1(Tab$time);
  
  # Survival step functions
  f0 = stepfun(x=K0$time,y=c(1,K0$surv));
  f1 = stepfun(x=K1$time,y=c(1,K1$surv));
  
  # Add survival to table
  Tab$surv0 = f0(Tab$time);
  Tab$surv1 = f1(Tab$time);
  
  # Add cumulative hazard
  Tab$cumhaz0 = cumsum(Tab$events0/Tab$nar0);
  Tab$cumhaz1 = cumsum(Tab$events1/Tab$nar1);
  
  # SE step functions
  f0 = stepfun(x=K0$time,y=c(0,K0$std.err));
  f1 = stepfun(x=K1$time,y=c(1,K1$std.err));
  
  # Add SEs to table
  Tab$se0 = f0(Tab$time);
  Tab$se1 = f1(Tab$time);
  
  # Hall-Wellner weights
  if(weights=="hw"){
    Tab$weight = 1/(1+Tab$se0^2+Tab$se1^2);
  } else {
    Tab$weight = 1;
  }
  
  # "Loop"
  aux = function(r){
    # Perturbations
    Z0 = sapply(Tab$events0,function(x){rnorm(n=1,mean=0,sd=sqrt(x))});
    Z1 = sapply(Tab$events1,function(x){rnorm(n=1,mean=0,sd=sqrt(x))});
    # Loop over doots
    V0 = (Z0/Tab$nar0);
    V1 = (Z1/Tab$nar1);
    # Cumsums
    Q0 = cumsum(V0);
    Q1 = cumsum(V1);
    # Scale
    R0 = Tab$surv0*Q0;
    R1 = Tab$surv1*Q1;
    # Scale
    R = -sqrt(n)*(Tab$weight)*(R0-R1);
    # Output
    return(max(abs(R)))
  }
  
  # Simulations
  Sups = sapply(seq(1:paths),aux);
  Ka = as.numeric(quantile(Sups,probs=1-sig));
  
  # Results
  Tab$D = (Tab$surv1-Tab$surv0);
  Tab$L = Tab$D-Ka/(Tab$weight*sqrt(n));
  Tab$U = Tab$D+Ka/(Tab$weight*sqrt(n));
  
  # Output
  Out = new(Class="band",Alpha=sig,Crit=Ka,Paths=paths,Samples=2,Table=Tab);
  return(Out);
}

