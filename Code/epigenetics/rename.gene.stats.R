#This function renames the strains in a gene.stats object
#from epigen.states.and.expr for plotting

rename.gene.stats <- function(gene.stats, col.table){
	state.mat.names <- colnames(gene.stats[[2]])
	colnames(gene.stats[[2]]) <- col.table[match(state.mat.names, col.table[,4]),2]

	expr.names <- names(gene.stats[[3]])
	names(gene.stats[[3]]) <- col.table[match(substr(expr.names, 1, 2), col.table[,4]),2]
	
	chrom.tree.names <- gene.stats[[5]]$tip.label
	gene.stats[[5]]$tip.label <- col.table[match(chrom.tree.names, col.table[,4]),2]
	
	return(gene.stats)
	}
	
