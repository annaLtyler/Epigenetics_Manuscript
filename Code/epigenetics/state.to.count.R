#This function converts a state matrix to a count
#matrix. States that are present/absent get a 1,
#States that are missing (NA) get a 0.

state.to.count <- function(state.mat){
    state.mat[which(is.finite(state.mat))] <- 1
    state.mat[which(!is.finite(state.mat))] <- 0
    return(state.mat)
} 