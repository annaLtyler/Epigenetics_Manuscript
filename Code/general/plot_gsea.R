plot_gsea <- function(ordered.genes, go.genes, go.name = NULL, min.size = 10, 
max.size = 600, nperm = 100, p = 0){

  gsea.results <- gsea(ordered.genes, go.genes, min.size = min.size, 
    max.size = max.size, nperm = nperm, p = p)
  
  gsea.enrich <- gsea.results[[1]]

  if(length(gsea.enrich) == 1){
    plot.text("Number of Genes Outside Range")
    invisible(NA)
  }

  consec.score <- consec.pairs(gsea.enrich)
  score.diffs <- consec.score[,2] - consec.score[,1]
  pos.locale <- which(score.diffs > 0)

  if(nperm > 0){
    par(mfrow = c(1,2))
  }
  plot(gsea.enrich, main = paste(go.name, "\n", length(pos.locale), "Genes"), 
  type = "l")
  abline(h = 0)
  if(length(pos.locale) > 0){
    points(x = pos.locale, y = rep(0, length(pos.locale)), col = "red", pch = "|")
  }

  result <- list("Enrichment" = gsea.enrich)

  if(nperm > 0){
    gsea.perm <- gsea.results[[2]]
    hist(gsea.perm, breaks = 25, main = "Permuted ES")
    abline(v = max(gsea.enrich), col = "red")
    result$"Perm.ES" <- gsea.perm
  }
  
  invisible(result)
}