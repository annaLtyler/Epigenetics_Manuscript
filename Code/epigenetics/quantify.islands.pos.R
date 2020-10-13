#This function takes a methyl matrix and quantifies 
#the methylation islands across strains.
#methyl.mat should have strains in rows and methylation 
#percentages in columns
#unlike quantify islands, this function uses position
#rather than changes in methylation percent to define
#islands.
#gene.frac.gap is the fraction of the length of the gene that
#consitutes a gap between methylation islands.
#If there are no islands at this level, which there usually are,
#the function then looks for islands using the 95th percentile
#of the gap size.

quantify.islands.pos <- function(methyl.mat, gene.frac.gap = 0.02, plot.bins = FALSE){

    if(length(methyl.mat) == 1){
        return(NA)
    }

    pos <- as.numeric(colnames(methyl.mat))
    #plot(pos)    
    consec.pos <- consec.pairs(pos)
    inter.pos <- consec.pos[,2] - consec.pos[,1]
    #head(cbind(consec.pos, inter.pos), 30)
    #hist(inter.pos, breaks = 200)
    gene.len <- max(pos) - min(pos)
    big.break <- round(gene.len*gene.frac.gap)
    transition.pts <- which(inter.pos > big.break)

    if(length(transition.pts) == 0){
        #use the 95th percentile of the break sizes
        big.break <- get.percentile(inter.pos, 95)
        transition.pts <- which(inter.pos > big.break)
    }

    #====================================================
    # internal functions
    #====================================================
    assign.bins <- function(transition.pts){
        if(length(transition.pts) == 0){
            return(rep(1, ncol(methyl.mat)))
        }
        island.bins <- rep(1,ncol(methyl.mat))        
        num.transitions <- length(transition.pts)

        for(i in 1:num.transitions){
            start.bin <- transition.pts[i]+1
            island.bins[start.bin:length(island.bins)] <- (i+1)
        }

        return(island.bins)
    }

   
    #====================================================
    
    island.bins <- assign.bins(transition.pts)
    names(island.bins) <- colnames(methyl.mat)
    num.islands <- max(island.bins)

    #plot(island.bins, pos)

    island.locale <- lapply(1:num.islands, function(x) which(island.bins == x))
    strain.islands <- lapply(island.locale, function(x) methyl.mat[,x,drop=FALSE])

    island.avg.methyl <- lapply(strain.islands, function(x) rowMeans(x, na.rm = TRUE))
        
    if(plot.bins){
        #quartz(width = 9, height = 5)
        plot.methyl.mat(methyl.mat, gene.name, xlim = NULL, island.bins)
    }

    result <- list("island.position" = island.bins, "island.avg.methyl" = island.avg.methyl)
    return(result)

}