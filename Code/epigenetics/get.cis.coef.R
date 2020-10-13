#This function gets the cis coefficient
#from a gene given its index in do.eqtl
#from epigen.states.and.expr
get.cis.coef <- function(do.eqtl, index){

  	idx.info <- do.eqtl[[index]]
	coef.list <- idx.info$DO.coefficients
	if(!is.null(coef.list)){
		coords <- as.numeric(unlist(lapply(strsplit(rownames(coef.list), "_"), function(x) x[2])))
  		gene.pos <- mean(c(idx.info[[1]]$Start.Mbp, idx.info[[1]]$End.Mbp)) * 1e6
		closest.idx <- which.min(abs(coords - gene.pos))
  		cis.coef <- idx.info$DO.coefficients[closest.idx,]
	  	# plot(coords, idx.info[[3]], type = "l")
  		# abline(v = gene.pos)
	  	return(cis.coef)	
	  	}else{
  		return(rep(0, ncol(do.eqtl[[1]]$DO.coefficients)))
	  	}
  	}
