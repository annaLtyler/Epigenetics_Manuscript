#This function gets the LOD score for each gene
#at the marker nearest its TSS


get.cis.LOD <- function(do.eqtl.data){

    info.locale <- which(names(do.eqtl.data) == "gene.info")
    gene.loc <- do.eqtl.data[[info.locale]][,"Start.Mbp"]*1e6

    marker.pos <- as.numeric(sapply(strsplit(rownames(do.eqtl.data$DO.coefficients), "_"), function(x) x[2]))
    nearest.pt <- get.nearest.pt(marker.pos, gene.loc)

    nearest.lod <- do.eqtl.data$LOD[nearest.pt]
    return(nearest.lod)
    
}