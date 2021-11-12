#This function bins an x vector into the specified number of
#bins and plots the y values in each bin of x.

plot_in_bins <- function(x, y, n.bins = 10, bin.type = c("even", "percentile"), main, xlab, ylab){

	if(missing(main)){main = ""}
	if(missing(xlab)){xlab = deparse(substitute(x))}
	if(missing(ylab)){ylab = deparse(substitute(y))}

    if(bin.type == "even"){
        x.bins <- segment.region(min(x, na.rm = TRUE), max(x, na.rm = TRUE)*1.01, n.bins, "ends")
        bin.idx <- lapply(1:(length(x.bins)-1), 
            function(a) intersect(which(x >= x.bins[a]), which(x <= x.bins[(a+1)])))
        binned.y <- lapply(bin.idx, function(a) y[a])
        binned.x <- lapply(bin.idx, function(a) x[a])
    }else{
        x.bins <- chunkV(sort(x), num.chunks)
        binned.y <- lapply(x.bins, function(a) y[match(names(a), names(y))])
        binned.x <- lapply(x.bins, function(a) x[match(names(a), names(x))])
    }

    names(binned.y) <- sapply(binned.x, function(x) signif(mean(x, na.rm = TRUE),2))
    boxplot(binned.y, xlab = xlab, ylab = ylab, main = main)

    results <- list("x" = binned.x, "y" = binned.y)
    invisible(results)
}