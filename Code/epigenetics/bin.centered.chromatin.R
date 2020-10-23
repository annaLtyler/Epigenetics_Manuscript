#This function bins a chromatin matrix into a fixed number
#of bins. It is meant to bin chromatin matrices whose coordinates
#have been scaled to the gene body. Upstream coordinates are < 0,
#downstream coordinates are > 1, and gene body coordinates are 
#between 0 and 1. To compare state presence/absence across multiple 
#genes in a spatially-aware way, we bin all genes into a fixed number
#of bins. This function looks across those bins and tallies whether
#each state is present or absent in those bins. 
#To calculate a proportion of each state in each bin, use tally.type == "proportion"
#We generate one presence/absence matrix or proportion matrix for each 
#state


bin.centered.chromatin <- function(centered.chrom.mat, coord.min = -1, coord.max = 2,
nbins = 500, nstates = 8, tally.type = c("present", "proportion")){

    tally.type = tally.type[1]

    n.strains <- ncol(centered.chrom.mat)
    coords <- as.numeric(rownames(centered.chrom.mat))

    bins <- segment.region(coord.min, coord.max, nbins, "ends")
    consec.bins <- consec.pairs(bins)

    binned.state.mats <- vector(mode = "list", length = nstates)
    names(binned.state.mats) <- 1:nstates
    for(i in 1:nstates){
        binned.mat <- matrix(NA, ncol = ncol(centered.chrom.mat), nrow = length(bins))
        rownames(binned.mat) <- bins
        colnames(binned.mat) <- colnames(centered.chrom.mat)
        binned.state.mats[[i]] <- binned.mat
    }

    for(i in 1:nrow(consec.bins)){
        bin.min <- consec.bins[i,1]
        bin.max <- consec.bins[i,2]
        incoords <- intersect(which(coords >= bin.min), which(coords <= bin.max))
        if(length(incoords) > 0){
            chrom.states <- centered.chrom.mat[incoords,,drop=FALSE]
            if(tally.type == "present"){
                state.count <- sapply(1:ncol(centered.chrom.mat), function(y) 
                    sapply(1:nstates, function(x) as.logical(length(which(chrom.states[y] == x)))))
            }else{
                state.count <- sapply(1:ncol(centered.chrom.mat), function(y) 
                    sapply(1:nstates, function(x) length(which(chrom.states[y] == x))/length(incoords)))
            }
            #rownames(state.count) <- 1:nstates
            #colnames(state.count) <- colnames(centered.chrom.mat)

            for(s in 1:nstates){
                binned.state.mats[[s]][i,] <- as.numeric(state.count[s,])
            }
        }
    }

    #use adjacent rows with values to fill in rows with all missing values.
    for(m in 1:length(binned.state.mats)){
        check.mat <- binned.state.mats[[m]]
        #pheatmap(check.mat, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE)
        not.na.row <- which(!is.na(rowSums(check.mat)))
        min.row <- min(not.na.row)
        max.row <- max(not.na.row)
        all.na.row <- which(is.na(rowSums(check.mat)))

        to.fill <- all.na.row[intersect(which(all.na.row >= min.row), which(all.na.row <= max.row))]
        if(length(to.fill) > 0){ #fill each NA row with the nearest filled row
            for(i in 1:length(to.fill)){
                #fill.with <- not.na.row[get.nearest.pt(not.na.row, to.fill[i])]
                fill.with <- not.na.row[max(which(not.na.row < to.fill[i]))]
                check.mat[to.fill[i],] <- check.mat[fill.with,]
            }
        }
        #pheatmap(check.mat, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE)
        
        binned.state.mats[[m]] <- check.mat

    }

    return(binned.state.mats)

}