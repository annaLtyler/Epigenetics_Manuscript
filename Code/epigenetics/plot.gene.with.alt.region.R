#This function plots stats for genes in a list and 
#puts the results in a pdf with the given label
  	
 plot.gene.with.alt.region <- function(qtl.ids, do.eqtl, rnaseq.gene.info, rna.seq, gene.name, pdf.label, bg.list = NULL, start.feature, end.feature, start.position, end.position, upstream.buffer = 1000, downstream.buffer = 1000){
	pdf(paste("Chromatin_and_expression_alt_region_", pdf.label, "_", start.position, "_", end.position, ".pdf", sep = ""))


	# for(sn in 1:length(gene.list)){
		# report.progress(sn, length(gene.list))
		# print(gene.list[sn])
		# gene.name <- gene.list[sn]
		id <- unique(rnaseq.gene.info[which(rnaseq.gene.info[,"external_gene_name"] == gene.name),1])
		if(length(id) > 0){
			gene.stats <- get.alt.gene.info(rnaseq.gene.info, rna.seq, id, start.position, end.position, upstream.buffer, downstream.buffer)
			gene.stats <- rename.gene.stats(gene.stats, col.table)
			plot.expr.and.states(gene.info = gene.stats[[1]], state.mat = gene.stats[[2]], gene.expr = gene.stats[[3]], snp.tree = gene.stats[[4]], chrom.tree = gene.stats[[5]])
			gene.start <- gene.stats[[1]][1,4]
			gene.end <- gene.stats[[1]][1,5]
			qtl.locale <- which(qtl.ids == id)	
			if(length(qtl.locale) > 0){
				plot.eqtl(coef.list = do.eqtl[[qtl.locale]][[2]], lod.score = do.eqtl[[qtl.locale]][[3]], label = gene.name, gene.coord = mean(c(gene.start, gene.end)))
				}
			}
		# }
	dev.off()
}	
	
