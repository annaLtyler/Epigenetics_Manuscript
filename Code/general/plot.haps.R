plot.haps <- function(hap.mat, ylab = "Individual"){

    require(qtl2)
    data(CCcolors)

    rowdist <- dist(hap.mat)
    rowclust <- hclust(rowdist)
    clust.mat <- hap.mat[rowclust$order,]

    col.mat <- matrix("lightgray", ncol = ncol(clust.mat), nrow = nrow(clust.mat))
    for(i in 1:ncol(clust.mat)){
        col.mat[which(clust.mat[,i] > 0.2),i] <- CCcolors[i]
    }
		
    plot(c(1, dim(col.mat)[2]), c(1, dim(col.mat)[1]), type = "n", axes = FALSE, 
    xlab = "", ylab = "", xlim = c(0.7, dim(col.mat)[2]+0.2), 
    ylim = c(0.7, dim(col.mat)[1]+0.2), bg = "transparent")

    rasterImage(col.mat, xleft = 0.5, ybottom = 0.5, xright = dim(col.mat)[2]+0.5, 
    ytop = dim(col.mat)[1]+0.5, interpolate = FALSE, bg = "transparent")

    par(xpd = TRUE)
    text(x = 1:ncol(col.mat), y = rep(nrow(col.mat)*1.02, ncol(col.mat)), 
    colnames(hap.mat))
    text(x = 0, y = nrow(col.mat)/2, labels = ylab, srt = 90)
    par(xpd = FALSE)

}