#This gene gets the ChromHMM state proportions of 
#a given gene
get.chrom.state.prop <- function(chrom.state.mats, num.states = 6, verbose = FALSE){
	
		get.prop <- function(state.mat, num.states){
			if(!is.null(nrow(state.mat))){
			state.prop.mat <- apply(state.mat, 2, function(x) unlist(lapply(1:num.states, function(y) length(which(x == y))/length(x))))
			rownames(state.prop.mat) <- 1:num.states
			return(state.prop.mat)
			}else{
				return(NA)
				}
			}
			
		if(verbose){
			all.state.prop <- lapply_pb(chrom.state.mats, function(x) get.prop(x, num.states))
			}else{
			all.state.prop <- lapply(chrom.state.mats, function(x) get.prop(x, num.states))	
			}
		return(all.state.prop)
	}

