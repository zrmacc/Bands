# Purpose: 2-sample KS Test

#' 2-Sample Kolmogorov Smirnov Test
#' 
#' @param time Observation time.
#' @param status Status, 1 for an event, 0 for censoring.
#' @param arm Treatment indicator, 1 for target group, 0 for reference.
#' @param paths Sample paths. 
#' @param weights Either equiprecision "ep" or Hall-Wellner "hw". 
#' 
#' @return List containing the KS statistic and the p-value. 
#' @export
#' 
#' @importFrom data.table data.table 
#' @importFrom methods new
#' @importFrom stats quantile rnorm stepfun
#' @importFrom survival Surv survfit 

survKS = function(time,status,arm,paths=1e3,weights="ep"){
  # Sample size
  n = length(time);
  # KM Curves
  K0 = survfit(Surv(time[arm==0],status[arm==0])~1);
  K1 = survfit(Surv(time[arm==1],status[arm==1])~1);
  
  # Lower endpoint
  L = max(min(K0$time[K0$surv<1]),min(K1$time[K1$surv<1]));
  # Upper endpoint
  U = min(max(K0$time[K0$surv>0]),max(K1$time[K1$surv>0]));
  
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
  
  # SE step functions
  f0 = stepfun(x=K0$time,y=c(0,K0$std.err));
  f1 = stepfun(x=K1$time,y=c(1,K1$std.err));
  
  # Add SEs to table
  Tab$se0 = f0(Tab$time);
  Tab$se1 = f1(Tab$time);
  
  # Hall-Wellner weights
  if(weights=="ep"){
    Tab$weight = 1;
  } else if(weights=="hw"){
    Tab$weight = 1/(1+Tab$se0^2+Tab$se1^2);
  }
  
  # Observed KS statistic
  Tab$D = (Tab$surv1-Tab$surv0);
  Dobs = max(Tab$weight*abs(Tab$D));
  
  # "Loop"
  aux = function(r){
    # Perturbations
    Z0 = sapply(Tab$events0,function(x){rnorm(n=1,mean=0,sd=sqrt(x))});
    Z1 = sapply(Tab$events1,function(x){rnorm(n=1,mean=0,sd=sqrt(x))});
    # Loop over doots
    V = -1*Tab$weight*(Tab$surv0*Z0/Tab$nar0-Tab$surv1*Z1/Tab$nar1);
    # Cumulative sum
    Q = cumsum(V);
    # Output
    return(max(abs(Q)))
  }
  
  # Simulations
  Sups = sapply(seq(1:(paths-1)),aux);
  
  # p-value
  p = (1+sum(Sups>=Dobs))/(1+length(Sups));
  
  # Report
  cat("2 Sample KS Test.\n");
  cat("Test statistic:",round(Dobs,digits=2),"\n");
  cat("p value:",signif(p,digits=3),"\n");
  cat("Sample paths:",paths,"\n");
  
  # Output
  Out = list("Statistic"=Dobs,"p"=p);
  
  return(Out);
}