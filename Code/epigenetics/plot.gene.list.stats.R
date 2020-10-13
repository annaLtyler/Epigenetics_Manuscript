#This function plots stats for genes in a list and 
#puts the results in a pdf with the given label
  	
 plot.gene.list.stats <- function(qtl.ids, do.eqtl, rnaseq.gene.info, rna.seq, gene.list, list.label, bg.list = NULL, start.feature, end.feature, upstream.buffer, downstream.buffer){
	# pdf(paste("Chromatin_and_expression_", list.label, "_", start.feature, "_", end.feature, "_", upstream.buffer, "_", downstream.buffer, ".pdf", sep = ""))
	if(!is.null(bg.list)){
		enrichment <- gprofiler(gene.list, organism = "mmusculus", custom_bg = bg.list, hier_filtering = "moderate")
		}else{
		enrichment <- gprofiler(gene.list, organism = "mmusculus", hier_filtering = "moderate")	
		}
	plot.enrichment(enrichment, 20, 0.7, list.label)
	for(sn in 1:length(gene.list)){
		# report.progress(sn, length(gene.list))
		print(gene.list[sn])
		gene.name <- gene.list[sn]
		id <- unique(rnaseq.gene.info[which(rnaseq.gene.info[,"external_gene_name"] == gene.name),1])
		if(length(id) > 0){
			gene.stats <- get.gene.info(rnaseq.gene.info, rna.seq, id, start.feature, end.feature, upstream.buffer, downstream.buffer)
			gene.stats <- rename.gene.stats(gene.stats, col.table)
			plot.expr.and.states(gene.info = gene.stats[[1]], state.mat = gene.stats[[2]], gene.expr = gene.stats[[3]], snp.tree = gene.stats[[4]], chrom.tree = gene.stats[[5]])
			gene.start <- gene.stats[[1]][1,4]
			gene.end <- gene.stats[[1]][1,5]
			qtl.locale <- which(qtl.ids == id)	
			if(length(qtl.locale) > 0){
				plot.eqtl(coef.list = do.eqtl[[qtl.locale]][[2]], lod.score = do.eqtl[[qtl.locale]][[3]], label = gene.name, gene.coord = mean(c(gene.start, gene.end)))
				}
			}
		}
	# dev.off()
}	
	
