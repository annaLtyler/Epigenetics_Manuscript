# calculate expanded chromtin array from haplotypes


calc.chrom <- function(gene.name, transcript.info, transcript.haplotypes, 
  chrom.states, strain.key, perm.order = NULL){
    
  gene.id <- transcript.info[which(transcript.info[,"external_gene_name"] == gene.name),"ensembl_gene_id"][1]

    if(length(gene.id) == 0){
        return(NA)
    }
    
  chrom.order <- match.order(colnames(transcript.haplotypes[[1]]),
    rownames(chrom.states[[1]]), strain.key)

  gene.idx <- which(names(transcript.haplotypes) == gene.id)
  chrom.idx <- which(names(chrom.states) == gene.id)
  
  if(length(gene.idx) == 0){return(NA)} #no expression data

  haps <- transcript.haplotypes[[gene.idx]]
  if(length(haps) == 1){return(NA)}
  if(length(chrom.states[[chrom.idx]]) == 1){return(NA)}

  ref.chroms <- chrom.states[[chrom.idx]][chrom.order,,]

  #if we want a permuted imputed chromatin state, shuffle 
  #the relationship between founder strain and chromatin
  #state vector
  if(!is.null(perm.order)){
    ref.chroms <- ref.chroms[perm.order,,]
    rownames(ref.chroms) <- rownames(chrom.states[[chrom.idx]][chrom.order,,])
  }

  chrom.array <- array(NA, dim = c(nrow(haps), ncol(ref.chroms), dim(ref.chroms)[3]))
  rownames(chrom.array) <- rownames(haps)
  colnames(chrom.array) <- paste0("State", colnames(ref.chroms))
  dimnames(chrom.array)[[3]] <- dimnames(ref.chroms)[[3]]
  
  for(p in 1:dim(ref.chroms)[3]){ #for each position
    chrom.mat <- ref.chroms[,,p]
    chrom.array[,,p] <- haps[,,1] %*% chrom.mat
    }
  
  return(chrom.array)
}