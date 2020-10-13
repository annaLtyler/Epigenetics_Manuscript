#This function plots a methylation table for a gene
#showing the gene body information

plot.methyl.table <- function(methyl.table, gene.name, gene.info.table, 
show.exons = TRUE, gene.lwd = 3){
	
	
	plot.methyl.mat(t(methyl.table), plot.label = gene.name)
	#imageWithText(t(methyl.table), show.text = FALSE, row.names = colnames(methyl.table), 
	#main = gene.name, use.pheatmap.colors = TRUE)
	
	gene.locale <- which(gene.info.table[,2] == gene.name)
	tss <- unique(gene.info.table[gene.locale,"start_position"])
	tes <- unique(gene.info.table[gene.locale,"end_position"])
	strand <- unique(gene.info.table[gene.locale,"strand"])
	
	gene.exon.table <- gene.info.table[gene.locale,]

	par(xpd = TRUE)
	if(tss != tes){
		if(strand == 1){
			arrows(tss, 0, tes, 0, lwd = gene.lwd, length = 0.1)
			}else{
			arrows(tes, 0, tss, 0, lwd = gene.lwd, length = 0.1)	
			}
		}else{
			points(tes, 0, cex = 2, pch = "*")	
			}


	if(show.exons){
	for(i in 1:length(gene.locale)){
		exon.start <- gene.exon.table[i,"exon_chrom_start"]
		exon.end <- gene.exon.table[i,"exon_chrom_end"]
		segments(exon.start, 0, exon.end, 0, lwd = gene.lwd, col = "yellow")
		}#end looping through exons
	}

	par(xpd = FALSE)
}


