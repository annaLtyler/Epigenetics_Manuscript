#This function plots stats for genes in a list and 
#puts the results in a pdf with the given label
  	
plot.one.gene.genetic <- function(gene.name, rnaseq.gene.info, 
genoprobs, map, pheno){

    gene.idx <- which(rnaseq.gene.info[,"external_gene_name"] == gene.name)
    gene.id <- unique(rnaseq.gene.info[gene.idx,"ensembl_gene_id"])
    gene.chr <- unique(rnaseq.gene.info[gene.idx,"chromosome_name"])
    gene.pos <- unique(rnaseq.gene.info[gene.idx,"start_position"])

    expr.locale <- which(colnames(pheno) == gene.id)

    chr.locale <- which(names(genoprobs) == gene.chr)

    lod <- scan1(genoprobs[,chr.locale], pheno[,expr.locale,drop=FALSE])
    

    nearest.pt <- get.nearest.pt(map[[chr.locale]], gene.pos/1e6)
    allele.coef <- scan1coef(genoprobs[,chr.locale], pheno[,expr.locale,drop=FALSE])
    
    par(mfrow = c(2,1))
    plot(lod, map = map, main = gene.name)
    abline(v = gene.pos/1e6)
    plot_coefCC(allele.coef, map = map)
    abline(v = gene.pos/1e6)

    gene.stats <- list("LOD" = lod, "coef" = allele.coef, 
    "nearest.lod" = lod[nearest.pt], "nearest.coef" = allele.coef[nearest.pt,])

    invisible(gene.stats)
}	
	
