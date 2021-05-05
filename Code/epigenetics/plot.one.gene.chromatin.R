#This function plots stats for genes in a list and 
#puts the results in a pdf with the given label
  	
plot.one.gene.chromatin <- function(do.eqtl, rnaseq.gene.info, rna.seq, all.bed, 
col.table, gene.name, T.or.C, start.feature, end.feature, upstream.buffer, 
downstream.buffer, get.snps = FALSE, separate.windows = FALSE, 
dim.rd = c("mds", "state.prop"), state = NULL, total.states = 9, 
show.state.numbers = FALSE,state.weights = NULL){


	qtl.ids <- names(do.eqtl)
	id <- unique(rnaseq.gene.info[which(rnaseq.gene.info[,"external_gene_name"] == gene.name),1])
	
	if(length(id) == 0){
		plot.new()
		plot.window(xlim = c(0, 1) ,ylim = c(0,1))
		text(0.5, 0.5, paste("There are no eQTL data for", gene.name))	
		return(NULL)
	}

	if(length(id) > 0){
		gene.stats <- get.gene.info(rnaseq.gene.info, rna.seq, all.bed, id, T.or.C,
		start.feature, end.feature, upstream.buffer, downstream.buffer, 
		get.snps = get.snps)
		
		if(is.null(gene.stats)){
			plot.text("No data for gene")
			return(NULL)
		}

		#gene.stats$control.expression <- lapply(gene.stats$control.expression, 
		#function(x) log2(x+1))
		
		plot.expr.and.states(gene.info = gene.stats$gene.info, state.mat = 
		gene.stats$state.matrix, gene.expr = gene.stats$control.expression, 
		snp.tree = gene.stats$snp.tree, chrom.tree = gene.stats$chromatin.state.tree, 
		dim.rd = dim.rd, state = state, num.states = total.states, 
		show.states = show.state.numbers, state.weights = state.weights, 
		col.table = col.table)
		
		
	}
invisible(gene.stats)
}	
	
