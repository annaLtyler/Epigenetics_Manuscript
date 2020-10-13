get.island.chrom <- function(selected.island, chrom.mat){
    
    if(length(selected.island) < 2){
        return(NA)
    }

    island.pos <- as.numeric(names(selected.island[[1]]))
    chrom.pos <- as.numeric(rownames(chrom.mat))
    overlapping.chrom <- intersect(which(chrom.pos >= min(island.pos)), which(chrom.pos <= max(island.pos)))
    island.chrom <- chrom.mat[overlapping.chrom,]
    return(island.chrom)

}