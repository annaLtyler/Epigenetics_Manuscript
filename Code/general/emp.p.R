#empirical p value calculation for our 
#small distributions

emp.p <- function(obs.val, null.dist, sig.fig = 3){
  if(is.na(obs.val)){return(NA)}
  as.or.more.extreme <- which(abs(round(null.dist, sig.fig)) >= abs(round(obs.val, sig.fig)))
  p <- length(as.or.more.extreme)/length(null.dist)
  return(p)
}
