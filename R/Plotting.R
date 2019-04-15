# Purpose: Plotting Confidence Bands
# Updated: 19/03/20

#' Master Plotting Function
#' 
#' @param B Object of class band.

plotBands = function(B){
  # Check samples
  m = B@Samples;
  q = NULL;
  if(m==1){
    q = plotOneSample(B);
  } else if(m==2){
    q = plotTwoSample(B);
  }
  return(q);
}

#' Plot 1-Sample Confidence Band
#' 
#' @param B Object of class band.
#' @import ggplot2

plotOneSample = function(B){
  # Data
  df = B@Table;
  # Plot
  q = ggplot(data=df) + theme_bw();
  q = q + geom_ribbon(aes(x=time,ymin=L,ymax=U),fill="#0073C2FF",alpha=0.5);
  q = q + geom_step(aes(x=time,y=surv),color="#0073C2FF");
  q = q + xlab("Time") + ylab("Survival");
  return(q);
}

#' Plot 2-Sample Confidence Band
#' 
#' @param B Object of class band.
#' @import ggplot2

plotTwoSample = function(B){
  # Data
  df = B@Table;
  # Plot
  q = ggplot(data=df) + theme_bw();
  q = q + geom_step(aes(x=time,y=D),color="#0073C2FF");
  q = q + geom_ribbon(aes(x=time,ymin=L,ymax=U),fill="#0073C2FF",alpha=0.5);
  q = q + geom_hline(yintercept=0,linetype="dashed",color="gray");
  q = q + xlab("Time") + ylab(expression(Survival~Difference~(S[1]-S[0])));
  return(q);
}