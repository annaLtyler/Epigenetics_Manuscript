# calculate expanded chromtin array from haplotypes


calc.chrom <- function(gene.name, transcript.info, transcript.haplotypes, 
chrom.states, strain.key){
    gene.id <- transcript.info[which(transcript.info[,"external_gene_name"] == gene.name),"ensembl_gene_id"][1]

    if(length(gene.id) == 0){
        return(NA)
    }
    
  chrom.order <- match.order(colnames(transcript.haplotypes[[1]]),
  rownames(chrom.states[[1]]), strain.key)

  gene.idx <- which(names(transcript.haplotypes) == gene.id)
  haps <- transcript.haplotypes[[gene.idx]]
  if(length(haps) == 1){return(NA)}
  ref.chroms <- chrom.states[[gene.idx]][chrom.order,,]

  chrom.array <- array(NA, dim = c(nrow(haps), ncol(ref.chroms), dim(ref.chroms)[3]))
  rownames(chrom.array) <- rownames(haps)
  colnames(chrom.array) <- paste0("State", colnames(ref.chroms))
  dimnames(chrom.array)[[3]] <- dimnames(ref.chroms)[[3]]
  
  for(p in 1:dim(ref.chroms)[3]){ #for each position
    chrom.mat <- ref.chroms[,,p]

    for(i in 1:dim(chrom.array)[1]){ #for each individual
      ind.hap <- haps[i,,1]
      #barplot(ind.hap)
      mult.mat <- apply(chrom.mat, 2, function(x) x*ind.hap)
      #pheatmap(mult.mat, cluster_rows = FALSE, cluster_cols = FALSE)
      ind.state <- colSums(mult.mat)
      #barplot(ind.state)
      #print(signif(ind.state, 2))
      chrom.array[i,,p] <- ind.state
    }
  }
  
  return(chrom.array)
}