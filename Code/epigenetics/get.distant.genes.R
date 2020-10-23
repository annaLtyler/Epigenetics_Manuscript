#This function gets lists of genes that are at least
#the specified distance from the nearest gene upstream
#and downstream. We do this by identifying any genes that
#are closer than the buffer and eliminating any genes that
#are within the specified distance of any other genes.

get.distant.genes <- function(gene.info, nearest.gene = 1000){
    chr <- unique(gene.info[,"chromosome_name"])

    distant.genes <- vector(mode = "list", length = length(chr))
    names(distant.genes) <- chr

    end.pairs <- pair.matrix(c(2,3), ordered = TRUE, self.pairs = TRUE)

    for(ch in 1:length(chr)){
        cat("Scanning chromosome", chr[ch], "...\n")
        chr.locale <- which(gene.info[,"chromosome_name"] == chr[ch])
        chr.genes <- as.matrix(unique(gene.info[chr.locale,c("external_gene_name", "start_position", "end_position")]))

        #get pairwise distances to other genes on the chromosome
        all.dist <- vector(mode = "list", length = nrow(end.pairs))
        for(d in 1:nrow(end.pairs)){
            start.coord <- end.pairs[d,1]
            end.coord <- end.pairs[d,2]
            dist.mat <- apply(chr.genes, 1, function(x) abs(as.numeric(x[start.coord]) - as.numeric(chr.genes[,end.coord])))
            diag(dist.mat) <- NA
            all.dist[[d]] <- dist.mat
        }

        near.genes <- lapply(all.dist, function(x) apply(x, 1, function(y) which(y < nearest.gene)))

        #check for length 0. This means there are no genes on the chromosome
        #within the distance specified.
        which.zero <- which(sapply(near.genes, length) == 0)
        if(length(which.zero) > 0){
            for(z in 1:length(which.zero)){
                near.genes[[which.zero[[z]]]] <- vector(mode = "list", length = nrow(chr.genes))
            }
        }        

        total.near <- lapply(1:nrow(chr.genes), function(x) Reduce("union", list(near.genes[[1]][[x]], near.genes[[2]][[x]], near.genes[[3]][[x]], near.genes[[4]][[x]])))

        num.near <- sapply(total.near, length)
        distant.genes[[ch]] <- chr.genes[which(num.near == 0),"external_gene_name"]
    }

    return(distant.genes)

}