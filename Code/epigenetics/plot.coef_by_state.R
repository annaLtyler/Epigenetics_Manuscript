#This function plots coefficient results not from
#scan1 as plot.coef does, but from a set of linear 
#models when we look at each state independently.
#used in both 2.1_chromaprobs_by_state.Rmd
#transcript.info is a table of information about transcripts.
#It is stored in RNASeq_gene_info.RData. This is required
#if you want to plot gene information.

plot.coef_by_state <- function(coef.results, gene.id = NULL, transcript.info = NULL, 
coef.type = c("chromatin", "genotype")){
  
  num.states <- nrow(coef.results)
  coef.type <- coef.type[1]
  if(coef.type == "genotype"){
    cols <- data(CCcolors)
    legend.labels <- LETTERS[1:8]
  }else{
    cols <- colors.from.values(1:num.states, use.pheatmap.colors = TRUE)
    legend.labels <- paste("State", 1:num.states)
  }
  ymin <- min(coef.results, na.rm = TRUE);ymax = max(coef.results, na.rm = TRUE)
  coords <- suppressWarnings(as.numeric(colnames(coef.results)))
  if(all(is.na(coords))){
    coords <- as.numeric(sapply(strsplit(colnames(coef.results), "_"), function(x) x[2]))
  }
  xmin <- min(coords); xmax <- max(coords)
  plot.new()
  plot.window(xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  for(s in 1:nrow(coef.results)){
    points(coords, coef.results[s,], col = cols[s], type = "l", lwd = 3)
  }
  axis(1);axis(2)
  mtext("Genomic Position", side = 1, line = 2.5)
  mtext("Coefficient", side = 2, line = 2.5)
  abline(h = 0, col = "gray", lty = 2)
  
  if(!is.null(gene.id)){
    gene.locale <- which(transcript.info[,"ensembl_gene_id"] == gene.id)[1]
    gene.name <- transcript.info[gene.locale,"external_gene_name"]
    gene.strand <- transcript.info[gene.locale,"strand"]
    mtext(gene.name, side = 3)
    gene.start <- transcript.info[gene.locale, "start_position"]    
    gene.end <- transcript.info[gene.locale, "end_position"]
    
    plot.new()
    plot.window(xlim = c(xmin, xmax), ylim = c(0, 1))
    if(gene.strand == 1){
      arrows(x0 = gene.start, x1 = gene.end, y0 = 0.5)
    }else{
      arrows(x0 = gene.end, x1 = gene.start, y0 = 0.5)
    }
    
  }
}