#This function plots the information in an island.obj
#from quantify.islands.pos for diagnosing issues.

plot.island.obj <- function(island.obj, add = FALSE){

    island.id <- island.obj$island.position
    pos <- as.numeric(names(island.id))
    u_islands <- unique(island.id)
    n.ind <- length(island.obj$island.avg.methyl[[1]])

    island.locale <- lapply(u_islands, function(x) which(island.id == x))
    
    if(!add){
    	    layout(matrix(c(1,2), nrow = 2), heights = c(1, 0.5))
	    par(mar = c(0,4,4,4))
    		plot.new()
	    plot.window(xlim = c(min(pos), max(pos)), ylim = c(0, n.ind))
	    }
    for(i in 1:length(u_islands)){
        draw.rectangle(min(pos[island.locale[[i]]]), max(pos[island.locale[[i]]]), 1, n.ind)
    }

	if(!add){
	    par(mar = c(4,4,0,4))
    		mean.vals <- sapply(island.obj[[2]], function(x) mean(x, na.rm = TRUE))
	    mean.pos <- sapply(island.locale, function(x) mean(as.numeric(names(x))))
    		plot.new()
	    plot.window(xlim = c(min(pos), max(pos)), ylim = c(0, 100))
    		segments(x0 = mean.pos, x1 = mean.pos, y0 = 0, y1 = mean.vals, lwd = 4)
	    axis(1);axis(2)
	    }

}