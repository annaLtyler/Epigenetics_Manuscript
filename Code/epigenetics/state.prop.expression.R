#This function correlates a state proportion with 
#inbred expression

state.prop.expression <- function(chrom.state.prop, group.mean.expr, strain.key, 
verbose = FALSE, plot.results = FALSE){
	
	num.states <- nrow(chrom.state.prop[[1]])
	null.result <- list("r" = rep(NA, num.states), "p" = rep(NA, num.states))
	strain.chrom.order <- match.order(strain.key[,1], colnames(chrom.state.prop[[1]]), strain.key)
	strain.expr.order <- match.order(strain.key[,1], names(group.mean.expr[[1]]), strain.key)

	one.cor <- function(gene.id){
		gene.id.chrom <- which(names(chrom.state.prop) == gene.id)

		#This case shouldn't really happen, if it does, 
		#we should probably rerun get.chrom.state
		if(length(gene.id.chrom) == 0){return(null.result)}

		gene.id.expr <- which(names(group.mean.expr) == gene.id)

		if(is.na(chrom.state.prop)[[gene.id.chrom]]){return(null.result)}
		
		gene.chrom <- chrom.state.prop[[gene.id.chrom]][,strain.chrom.order]
		gene.expr <- scale(group.mean.expr[[gene.id.expr]][strain.expr.order])
		#gene.expr <- group.mean.expr[[gene.id.expr]][strain.expr.order]

		if(all(is.na(gene.expr))){return(null.result)}

		if(plot.results){
			annot_row <- data.frame(gene.expr)
			colnames(annot_row) <- "Expression"
			rownames(annot_row) <- colnames(gene.chrom)
			pheatmap(t(gene.chrom), annotation_row = annot_row)
			#pheatmap(t(gene.chrom[,order(gene.expr)]), cluster_rows = FALSE, annotation_row = annot_row)
			par(mfrow = c(3,3))
			apply(gene.chrom, 1, function(x) plot(x, gene.expr))
		}

		all.cor <- apply(gene.chrom, 1, function(x) suppressWarnings(cor.test(gene.expr,x)))
		
		all.r  <- sapply(all.cor, function(x) x$estimate)		
		all.p <- sapply(all.cor, function(x) x$p.value)

		return(list("r" = all.r, "p" = all.p))
			
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