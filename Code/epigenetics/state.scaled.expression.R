#This function correlates a scaled state matrix from get.scaled.chrom.mats with
#inbred expression

state.scaled.expression <- function(scaled.chrom, group.mean.expr, strain.key, verbose = FALSE){
	
	null.result <- list("r" = NA, "p" = NA)
	strain.chrom.order <- match.order(strain.key[,1], rownames(scaled.chrom[[1]]), strain.key)
	strain.expr.order <- match.order(strain.key[,1], rownames(group.mean.expr[[1]]), strain.key)

	one.cor <- function(gene.id){
		gene.id.chrom <- which(names(scaled.chrom) == gene.id)

		#This case shouldn't really happen, if it does, 
		#we should probably rerun get.chrom.state
		if(length(gene.id.chrom) == 0){return(null.result)}

		gene.id.expr <- which(names(group.mean.expr) == gene.id)

		if(is.na(scaled.chrom)[[gene.id.chrom]]){return(null.result)}
		
		gene.chrom <- scaled.chrom[[gene.id.chrom]][strain.chrom.order,]
		
		if(length(gene.chrom) == 0){return(null.result)}
		
		gene.expr <- scale(group.mean.expr[[gene.id.expr]][strain.expr.order])

		if(all(is.na(gene.expr))){return(null.result)}

		all.cor <- cor.test(gene.expr,gene.chrom)
		all.r  <- all.cor$estimate
		all.p <- all.cor$p.value
		return(c("r" = all.r, "p" = all.p))
			
		}
	
	#for debugging
	#test.result <- for(i in 1:length(group.mean.expr)){one.cor(names(group.mean.expr)[i])}
	#gene.id <- names(group.mean.expr)[i]

	if(verbose){
		results <- lapply_pb(names(group.mean.expr), one.cor)
		}else{
		results <- lapply(names(group.mean.expr), one.cor)	
		}

	names(results) <- names(group.mean.expr)
	return(results)
	
}