#This function can be used with Reduce to find the running sum of 
#squares of a list of matrices.

ssq.mat <- function(A, B){
    AB <- ifelse(is.na(A),ifelse(is.na(B),NA,B), ifelse(is.na(B), A, A+(B^2))) 
    return(AB)
}
