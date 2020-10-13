chromState2Prop <- function(chrom.states, num.states, num.strains){
	
	if(is.null(chrom.states) || is.na(chrom.states)){
		return(matrix(NA, nrow = num.states, ncol = num.strains))
		}
	state.prop.mat <- apply(chrom.states, 2, function(x) unlist(lapply(1:num.states, function(y) length(which(x == y))/length(x))))
	rownames(state.prop.mat) <- 1:num.states
	return(state.prop.mat)
	
}