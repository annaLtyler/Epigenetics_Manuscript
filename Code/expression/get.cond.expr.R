#This function retrieves control expression 
#for a given gene

get.cond.expr <- function(rnaseq.gene.info, rna.seq, gene.id, T.or.C = "C",
return.mean = TRUE){
	
	gene.locale <- which(rnaseq.gene.info[,1] == gene.id)
	gene.expr.locale <- which(rownames(rna.seq) == gene.id)
	
	if(length(gene.locale) > 0){
		#pull out expression values for controls
		gene.expr <- rna.seq[gene.expr.locale,]	
		labels <- substr(names(gene.expr), 1, 3)
		u_labels <- sort(unique(labels))
		binned.expr <- lapply(u_labels, function(x) as.vector(gene.expr[which(labels == x)]))
		names(binned.expr) <- u_labels
		c.locale <- which(substr(names(binned.expr), 3, 3) == T.or.C)
		c.expr <- binned.expr[c.locale]
		ind.expr <- unlist(c.expr)
		expr.table <- Reduce("cbind", c.expr)
		colnames(expr.table) <- names(c.expr)
		if(return.mean){
			groups <- substr(names(ind.expr), 1, 2)
            u_groups <- unique(groups)
            strain.means <- unlist(lapply(u_groups, 
            function(x) mean(ind.expr[which(groups == x)])))
            names(strain.means) <- u_groups
		}else{
			strain.means <- expr.table
		}
	}else{
		strain.means <- NA
	}
return(strain.means)
}
