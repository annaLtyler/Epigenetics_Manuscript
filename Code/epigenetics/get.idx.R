
#This function finds the nearest index
#of a vector V to a given value val
get.idx <- function(V, val){
	val <- as.numeric(val)
	V <- as.numeric(V)

	if(val > max(V)){
		return(Inf)
		}
	if(val < min(V)){
		return(-Inf)
		}

	diffs <- V - val
	nearest.locale <- which(abs(diffs) == min(abs(diffs)))
	return(nearest.locale[1])
	}

