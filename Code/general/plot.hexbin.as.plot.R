#This plots a hexbin plot using a regular plotting function so 
#we can get multiple plots on one page. I'm still not sure how
#do to it when plotting an actual hexbin object.

plot.hexbin.as.plot <- function(x, y, xlab = "x", ylab = "y", main = "", 
    min.cex = 1, max.cex = 3, n.bins = 10, count.scale.fun = NULL,
    legend.pos = "topright", round.legend = 1000, use.pheatmap.colors = TRUE, 
    col.scale = "blue", grad.dir = "high", col.fun = c("linear", "exponential"), 
    exp.steepness = 1){

    if(missing(main)){main = ""}
	if(missing(xlab)){xlab = deparse(substitute(x))}
	if(missing(ylab)){ylab = deparse(substitute(y))}
    col.fun <- col.fun[1]

    hb <- hexbin(x, y, xlab = xlab, ylab = ylab)
    #plot(hb)
    xy <- hcell2xy(hb)
    pt.cex = scale.between.vals(hb@count, target.min = min.cex, target.max = max.cex)
    if(!is.null(count.scale.fun)){
        scale.fun <- match.fun(count.scale.fun)
        cols <- colors.from.values(scale.fun(hb@count), use.pheatmap.colors = use.pheatmap.colors,
        col.scale = col.scale, grad.dir = grad.dir, color.fun = col.fun, exp.steepness = exp.steepness)
    }else{
    cols <- colors.from.values(hb@count, use.pheatmap.colors = use.pheatmap.colors,
        col.scale = col.scale, grad.dir = grad.dir, color.fun = col.fun, exp.steepness = exp.steepness)
    }
    
    plot(xy, xlab = xlab, 
        ylab = ylab, pch = 16,
        main = main, col = cols, 
        cex = pt.cex)
    
    bins <- sort(unique(hb@count))
    binned.bins <- round(segment.region(min(bins), max(bins), n.bins, "ends"))
    rounded.bins <- unique(round(binned.bins/round.legend) * round.legend)
    rounded.bins[which(rounded.bins == 0)] <- 1

    bin.cols <- colors.from.values(rounded.bins, use.pheatmap.colors = use.pheatmap.colors,
    col.scale = col.scale, grad.dir = grad.dir, color.fun = col.fun, exp.steepness = exp.steepness)

    bin.cex <- scale.between.vals(rounded.bins, target.min = min.cex, target.max = max.cex)
    legend(legend.pos, legend = rounded.bins, col = bin.cols, pch = 16)

}