#This function takes the output of an EdgeR analysis (qlf object)
#and looks for functional enrichment in up and down regulated genes.
#It returns the up and down enrichments in a list for further processing.
#For example, plot.enrichment can be used to plot summary tables in 
#an easy-to-read format.

expr.fun.enrich <- function(qlf.obj, pval = 0.05, lFC = 2){

	 
	subtable <- qlf.obj$table
	if(!is.null(lFC)){
	  	subtable <- subtable[which(abs(subtable$logFC) >= lFC),]
	  	}
	if(!is.null(pval)){
		subtable <- subtable[which(subtable$PValue <= pval),]
	  	}
	
	downreg <- which(subtable[,"logFC"] < 0)
	upreg <- which(subtable[,"logFC"] > 0)
	
	get.enrichment <- function(gene.ids){
		gene.table <- gconvert(gene.ids, organism = "mmusculus")
		gene.names <- as.vector(gene.table[,"name"])
		enrichment <- gprofiler(gene.names, organism = "mmusculus")
		return(enrichment)
		}
	
	down.enrich <- get.enrichment(gene.ids = rownames(subtable)[downreg])
	up.enrich <- get.enrichment(gene.ids = rownames(subtable)[upreg])
	
	# plot.enrichment(down.enrich, text.size = 0.8, plot.label = "Downregulated with Treatment")
	# plot.enrichment(up.enrich, text.size = 0.8, plot.label = "Upregulated with Treatment")
	
	results <- list("downregulated" = down.enrich, "upregulated" = up.enrich)	
	return(results)
}