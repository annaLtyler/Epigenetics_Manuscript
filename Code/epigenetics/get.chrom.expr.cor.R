#This function starts with a gene id and 
#returns chromatin by expression correlation


get.chrom.expr.cor <- function(rnaseq.gene.info, rna.seq, gene.id, start.feature, end.feature, upstream.buffer, downstream.buffer){
	
	gene.locale <- which(rnaseq.gene.info[,1] == gene.id)
	gene.expr.locale <- which(rownames(rna.seq) == gene.id)
	
	if(length(gene.locale) > 0){
		gene.info <- rnaseq.gene.info[gene.locale,,drop=FALSE]
	
		#use larger buffers then trim to make sure we get states
		strain.states <- lapply(all.bed, function(x) get.states(gene.info, x, upstream.buffer = 50000, downstream.buffer = 50000, start.feature = start.feature, end.feature = end.feature))
		names(strain.states) <- strain.names
			
		if(sum(unlist(lapply(strain.states, nrow))) > 0){

			#pull out expression values for controls
			gene.expr <- rna.seq[gene.expr.locale,]	
			labels <- substr(names(gene.expr), 1, 3)
			u_labels <- sort(unique(labels))
			binned.expr <- lapply(u_labels, function(x) as.vector(gene.expr[which(labels == x)]))
			names(binned.expr) <- u_labels
			c.locale <- which(substr(names(binned.expr), 3, 3) == "C")
			c.expr <- binned.expr[c.locale]
			med.order <- order(unlist(lapply(binned.expr[c.locale], median)))

			
			#line up the states for the different strains and trim to requested region
			gene.mat <- align.states(strain.states)
			gene.mat <- trim.state.mat(gene.mat, gene.info, upstream.buffer, downstream.buffer, start.feature, end.feature)
	
			#get distance matrix for chromatin states and SNPs
			d.chrom.mat <- get.hamming.mat(gene.mat)
			#test whether chromatin state is related to expression
			expr.state.comp <- compare.expr.chrom(d.chrom.mat, c.expr)
			}else{
			expr.state.comp <- NA
			}
		}else{
			expr.state.comp <- NA
			}
	return(expr.state.comp)
}
