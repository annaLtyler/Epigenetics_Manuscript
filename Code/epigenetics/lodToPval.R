#This function converts LOD scores to p values

lodToPval <- function(x){
  pchisq(x*(2*log(10)),df=1,lower.tail=FALSE)/2
  }	
