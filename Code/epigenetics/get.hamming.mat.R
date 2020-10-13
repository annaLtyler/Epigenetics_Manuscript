#This function returns a Hamming distance matrix
#given a matrix of ChromHMM states for different 
#strains

get.hamming.mat <- function(gene.mat, na.rm = TRUE){
	#take out any rows with NAs. This isn't a great
	#solution. In the future, make sure 
	if(class(gene.mat)[1] != "matrix"){
		temp.mat <- matrix(gene.mat, nrow = 1)
		colnames(temp.mat) <- names(gene.mat)
		gene.mat <- temp.mat
		}

	na.locale <- which(is.na(gene.mat))
	
	if(length(na.locale) > 0){
		if(na.rm){
			not.na.rows <- which(!is.na(rowSums(gene.mat)))
			gene.mat <- gene.mat[not.na.rows,]
		}else{
			gene.mat[na.locale] <- 0
		}
	}

	
	if(is.null(ncol(gene.mat)) || nrow(gene.mat) == 0){
		return(NULL)
		}
	
	strain.pairs <- pair.matrix(1:ncol(gene.mat))	
	d.mat <- matrix(0, nrow = ncol(gene.mat), ncol = ncol(gene.mat))
	for(j in 1:nrow(strain.pairs)){
		v.dist <- hamming.distance(gene.mat[,strain.pairs[j,1]], gene.mat[,strain.pairs[j,2]])
		d.mat[strain.pairs[j,1], strain.pairs[j,2]] <- d.mat[strain.pairs[j,2], strain.pairs[j,1]] <- v.dist
		}
	return(d.mat)
	}
