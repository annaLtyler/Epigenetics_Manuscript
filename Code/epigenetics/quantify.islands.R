#This function takes a methyl matrix and quantifies 
#the methylation islands across strains.
#methyl.mat should have strains in rows and methylation 
#percentages in columns
#min.representation sets how many of the strains need to 
#be represented within a given island. If there are fewer
#than this number, the islands are merged with neighboring 
#islands


quantify.islands <- function(methyl.mat, transition.percent = 50, min.representation = 3,
plot.peaks = FALSE){

    #pheatmap(methyl.mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE)

    #avg.methyl <- apply(methyl.mat, 2, function(x) mean(x, na.rm = TRUE))
    min.methyl <- apply(methyl.mat, 2, function(x) min(x, na.rm = TRUE))
    max.methyl <- apply(methyl.mat, 2, function(x) max(x, na.rm = TRUE))
    methyl.var <- apply(methyl.mat, 2, function(x) var(x, na.rm = TRUE))
    methyl.var[which(is.na(methyl.var))] <- 0 

   if(plot.peaks){
        quartz(width = 11, height = 6)
        #par(mfrow = c(2,1))
        plot.new();plot.window(xlim = c(1,length(min.methyl)), ylim = c(0, 100))
        #points(avg.methyl, type = "l")
        points(min.methyl, type = "l", col = "red")
        points(max.methyl, type = "l", col = "blue")
        abline(h = transition.percent)
        axis(1);axis(2)
        #plot(methyl.var, type = "l")
   }
    
    min.slopes <- consec.slope(1:length(min.methyl), min.methyl)
    max.slopes <- consec.slope(1:length(max.methyl), max.methyl)
    var.slopes <- consec.slope(1:length(methyl.var), methyl.var)
    #plot(min.slopes, max.slopes)

    transition.pts <- intersect(which(abs(max.slopes) >= transition.percent), 
    which(abs(min.slopes) >= transition.percent))+1

    #percent.var <- (methyl.var/max(abs(methyl.var)))*100
    #transition.pts <- which(percent.var > transition.percent)

    if(length(transition.pts) == 0){
        transition.pts <- sort(union(which(abs(max.slopes) >= transition.percent), 
        which(abs(min.slopes) >= transition.percent)))+1
    }

    if(length(transition.pts) == 0){
        stop("There are no obvious places to bin the methylation data. Try reducing the transition percent.")
    }

    #====================================================
    # internal functions
    #====================================================
    assign.bins <- function(transition.pts){
        if(length(transition.pts) == 0){
            return(rep(1, length(min.methyl)))
        }
        island.bins <- rep(1,length(min.methyl))
        bin.num <- 2
        num.transitions <- max(c(1, length(transition.pts)-1))
        for(i in 1:num.transitions){
            start.bin <- transition.pts[i]
            if(num.transitions == 1){
                end.bin <- length(island.bins)   
            }else{
                end.bin <- transition.pts[i+1]
            }
            island.bins[start.bin:end.bin] <- bin.num
            bin.num <- bin.num + 1
        }

        if(num.transitions > 1){
            island.bins[end.bin:length(island.bins)] <- bin.num
        }
        return(island.bins)
    }

    remove.transitions <- function(transition.pts, low.island.idx){
        last.island <- any(low.island.idx > length(transition.pts))
        transition.pts <- transition.pts[-low.island.idx]
        if(last.island){
            transition.pts <- transition.pts[-length(transition.pts)]
        }
        return(transition.pts)
    }
    #====================================================
    
    
    island.bins <- assign.bins(transition.pts)
    num.islands <- max(island.bins)

    island.locale <- lapply(1:num.islands, function(x) which(island.bins == x))
    strain.islands <- lapply(island.locale, function(x) methyl.mat[,x,drop=FALSE])

    #The following line not only finds islands with low representation
    #But also islands with no methylation across multiple strains.
    num.strain.per.island  <- sapply(strain.islands, function(x) length(which(rowSums(x, na.rm = TRUE) > 0)))
    low.island <- which(num.strain.per.island < min.representation)

    #Islands with low representation or no methylation are merged with
    #neighboring islands
    #if there are small islands, remove these 
    while(length(low.island) > 0){
        transition.pts <- remove.transitions(transition.pts, low.island)
        island.bins <- assign.bins(transition.pts)
        
        num.islands <- max(island.bins)
        island.locale <- lapply(1:num.islands, function(x) which(island.bins == x))
        strain.islands <- lapply(island.locale, function(x) methyl.mat[,x,drop=FALSE])
        
        num.strain.per.island  <- sapply(strain.islands, function(x) length(which(rowSums(x, na.rm = TRUE) > 0)))
        low.island <- which(num.strain.per.island < min.representation)
    }

    island.avg.methyl <- lapply(strain.islands, function(x) rowMeans(x, na.rm = TRUE))
        
    result <- list("island.position" = island.bins, "island.avg.methyl" = island.avg.methyl)
    return(result)

}