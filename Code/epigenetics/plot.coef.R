#This function plots coefficient results from scan1
#in both 2.1_chromaprobs.Rmd and 2.2_chromaprobs_analysis.Rmd
#transcript.info is a table of information about transcripts.
#It is stored in RNASeq_gene_info.RData. This is required
#if you want to plot gene information.

plot.coef <- function(scan1.results, gene.id = NULL, transcript.info = NULL, 
coef.type = c("chromatin", "genotype")){
  
  require(qtl2)
  require(pheatmap)

  coef.type <- coef.type[1]
  if(coef.type == "genotype"){
    cols <- data(CCcolors)
    legend.labels <- LETTERS[1:8]
  }else{
    cols <- colors.from.values(1:8, use.pheatmap.colors = TRUE)
    legend.labels <- paste("State", 1:8)
  }
  just.coef <- scan1.results[,1:8]
  ymin <- min(just.coef, na.rm = TRUE);ymax = max(just.coef, na.rm = TRUE)
  coords <- suppressWarnings(as.numeric(rownames(just.coef)))
  if(all(is.na(coords))){
    coords <- as.numeric(sapply(strsplit(rownames(just.coef), "_"), function(x) x[2]))
  }
  xmin <- min(coords); xmax <- max(coords)
  plot.new()
  plot.window(xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  for(s in 1:ncol(just.coef)){
    points(coords, just.coef[,s], col = cols[s], type = "l", lwd = 3)
  }
  axis(1);axis(2)
  mtext("Genomic Position", side = 1, line = 2.5)
  mtext("Coefficient", side = 2, line = 2.5)
  if(!is.null(gene.id)){
    gene.locale <- which(transcript.info[,"ensembl_gene_id"] == gene.id)[1]
    gene.name <- transcript.info[gene.locale,"external_gene_name"]
    mtext(gene.name, side = 3)
    gene.start <- transcript.info[gene.locale, "start_position"]    
    gene.end <- transcript.info[gene.locale, "end_position"]    
    abline(v = c(gene.start, gene.end))
  }
}