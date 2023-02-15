#from https://stats.stackexchange.com/questions/6940/calculating-p-value-from-an-arbitrary-distribution

emp.p <- function(obs.val, null.dist, sig.fig = 3){
  if(is.na(obs.val)){return(NA)}
  
  dF <- function(null.dist){
    dnorm(null.dist)
  }
  
  pF <- function(q){
    integrate(dF,-Inf,q)$value
  }

  p <- pF(obs.val)

  return(p)
}
