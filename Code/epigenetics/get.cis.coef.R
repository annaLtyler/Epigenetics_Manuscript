#This function gets the cis coefficient
#from a gene given its index in do.eqtl
#from epigen.states.and.expr
get.cis.coef <- function(do.eqtl, mrna.info, index){

	if(is.null(do.eqtl[[index]])){
		return(rep(NA, 3))
	}

  	transcript.id <- names(do.eqtl)[index]
	nearest_marker <- unlist(mrna.info[index,"nearest.marker.id"])
	if(!is.null(nearest_marker) && !is.na(nearest_marker)){
		marker_idx <- which(rownames(do.eqtl[[index]]) == nearest_marker)
		cis.coef <- do.eqtl[[index]][marker_idx]
		result <- c(transcript.id, nearest_marker, cis.coef)
	}else{
		result <- c(transcript.id, nearest_marker, NA)
	}
	return(result)
  	}
