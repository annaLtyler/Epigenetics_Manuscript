#returns the variance explained of a random forest

rf.var <- function(rf.obj){
  var.exp <- tail(rf.obj$rsq, 1)*100
  if(length(var.exp) == 0){
    var.exp <- NA
  }
  return(var.exp)
}