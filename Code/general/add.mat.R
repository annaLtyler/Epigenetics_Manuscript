#This function can be used with Reduce to add a list of matrices
#together while ignoring NA entries.

add.mat <- function(A, B){
    AB <- ifelse(is.na(A),ifelse(is.na(B),NA,B), ifelse(is.na(B), A, A+B)) 
    return(AB)
}
