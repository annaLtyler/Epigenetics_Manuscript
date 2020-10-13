#This function creates a chromatin state matrix based on the
#methylation matrix. The final matched.chrom.mat has the same
#dimensions and coordinates as the methyl.mat, but instead of
#percentages, has the chromatin state at each position. 
#when positions are on chromatin state borders, the first of the
#two states is used.
#This function also produces a set of matrices corresponding to each
#chromatin state (1 through 8). Each matrix is the same dimensions
#as the original methyl.mat, and stores the methylation percentage
#at each position occupied by the given state. If there are no instances
#of state 2 in the matrix, the matrix for state 2 will be all NAs.


methylation.vs.state <- function(chrom.mat, methyl.mat, gene.name = NULL, gene.info = NULL,
plot.results = FALSE){

    if(length(methyl.mat) == 1 || length(chrom.mat) == 1){
        return(NA)
    }

    if(!is.null(gene.name)){
        gene.locale <- which(gene.info[,"external_gene_name"] == gene.name)
        gene.strand <- gene.info[gene.locale[1], "strand"]
    }

    chrom.pos <- as.numeric(rownames(chrom.mat))
    methyl.pos <- as.numeric(rownames(methyl.mat))

    #Chromatin states are in 200bp bins. Methylation
    #is at the bp level. For each methylation position,
    #calculate the chromatin state at that position

    match.chrom.state <- function(methyl.i){
        less.than <- which(chrom.pos < methyl.i)
        if(length(less.than) == 0){
            chrom.state <- rep(NA, ncol(methyl.mat))
        }
        greater.than <- which(chrom.pos > methyl.i)
        if(length(greater.than) == 0){
            chrom.state <- rep(NA, ncol(methyl.mat))
        }
        if(length(greater.than) > 0 && length(less.than) > 0){
            chrom.state <- chrom.mat[max(less.than),]
        }
        return(chrom.state)
    }

    matched.chrom <- sapply(methyl.pos, match.chrom.state)
    colnames(matched.chrom) <- methyl.pos
    rownames(matched.chrom) <- colnames(chrom.mat)

    #also make a methylation percent by position and state for each strain
    methyl.by.state.strain <- vector(mode = "list", length = 8)
    names(methyl.by.state.strain) <- 1:8
    for(i in 1:8){ #i is for the state number
        methyl.by.state.mat <- matrix(NA, nrow = ncol(methyl.mat), ncol = nrow(methyl.mat))
        colnames(methyl.by.state.mat) <- rownames(methyl.mat)
        rownames(methyl.by.state.mat) <- colnames(chrom.mat)
        for(j in 1:ncol(methyl.mat)){ #j is for strain
            state.locale <- which(matched.chrom[j,] == i)
            if(length(state.locale) > 0){
                state.methyl <- methyl.mat[state.locale,j]
                methyl.by.state.mat[j,state.locale] <- state.methyl
            }
        }
        methyl.by.state.strain[[i]] <- methyl.by.state.mat
    }

    #pheatmap(methyl.by.state.strain[[7]], cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE)


    if(plot.results){
        if(gene.strand == 1){
            pos.col <- colors.from.values(methyl.pos, use.pheatmap.colors = TRUE)
        }else{
            rev.pos <- (methyl.pos - max(methyl.pos))*-1
            pos.col <- colors.from.values(rev.pos, use.pheatmap.colors = TRUE)
        }
        par(mfrow = c(3,3))
        for(i in 1:ncol(methyl.mat)){
            plot(matched.chrom[i,], methyl.mat[,i],
            main = paste("Methylation By Chromatin State\n", gene.name, colnames(chrom.mat)[i]),
            xlab = "Chromatin State", ylab = "Methylation", col = pos.col)
        }
    }

    invisible(list(matched.chrom, methyl.by.state.strain))

}