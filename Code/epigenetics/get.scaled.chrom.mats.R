#This function performs multi-dimensional scaling on chromatin state matrices
#chrom.state.mats is the list of chromatin state matrices for each transcript
#d is the number of dimensions

get.scaled.chrom.mats <- function(chrom.state.mats, d = 1){
	
	all.hamming <- lapply_pb(chrom.state.mats, function(x) if(length(x) > 1){get.hamming.mat(x)}else{NA})	
	names(all.hamming) <- names(chrom.state.mats)		


	get.mds <- function(d.mat, d){
		if(length(d.mat) > 1){
			result <- cmdscale(d.mat, d)
			rownames(result) <- colnames(chrom.state.mats[[1]])
			return(result)
		}
		return(NA)
	}

	scaled.v <- lapply(all.hamming, function(x) get.mds(x, d))

	return(scaled.v)	
	
}