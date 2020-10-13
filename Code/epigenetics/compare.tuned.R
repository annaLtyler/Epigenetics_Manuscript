#This function takes in two objects returned from tune.maxnodes
#and plots them on the same axes to compare.

compare.tuned <- function(tuned.list, plot.label = "", plot.results = TRUE, legend.x = 0.17){

    non.null <- which(sapply(tuned.list, length) > 0)
    tuned.list <- tuned.list[non.null]

    xmin <- 0; xmax <- nrow(tuned.list[[1]])
    ymin1 <- min(unlist(sapply(tuned.list, function(x) x[,1])))
    ymax1 <- max(unlist(sapply(tuned.list, function(x) x[,1])))
    ymin2 <- min(unlist(sapply(tuned.list, function(x) x[,2])))
    ymax2 <- max(unlist(sapply(tuned.list, function(x) x[,2])))

    min.error <- sapply(tuned.list, function(x) which.min(x[,1]))
    max.exp <- sapply(tuned.list, function(x) which.max(x[,2]))

    if(plot.results){
        error.cols <- brewer.pal(9, "Set1")
        exp.cols <- sapply(error.cols, lighten)

        par(mar = c(5,5,5,15))
        plot.new()
        plot.window(xlim = c(xmin, xmax), ylim = c(ymin1, ymax1))
        
        for(i in 1:length(tuned.list)){
            points(1:nrow(tuned.list[[i]]), tuned.list[[i]][,1], 
            type = "l", col = error.cols[i], lwd = 2)
        }
    
        par(new = TRUE)
        plot.new()
        plot.window(xlim = c(xmin, xmax), ylim = c(ymin2, ymax2))
        for(i in 1:length(tuned.list)){
            points(1:nrow(tuned.list[[i]]), tuned.list[[i]][,2], 
            type = "l", col = exp.cols[i], lwd = 2)
        }
   
        axis(side = 2)
        mtext(side = 2, line = 2.5, "MSE")
    
        axis(side = 4)
        mtext(side = 4, line = 2.5, "Variance Explained")
    
        par(xpd = TRUE)
        plot.dim <- par("usr")
        plot.height <- plot.dim[4] - plot.dim[3]
        plot.width <- plot.dim[2] - plot.dim[1]

        legend(x = (plot.dim[2]+(plot.width*legend.x)), y = plot.dim[4], 
        col = error.cols, lty = 1, legend = names(tuned.list), 
        horiz = FALSE, lwd = 2,title = "MSE", cex = 0.8)

        legend(x = (plot.dim[2]+(plot.width*legend.x)), y = plot.dim[4]/2, 
        col = exp.cols, lty = 1, legend = names(tuned.list), 
        horiz = FALSE, lwd = 2,title = "Variance Explained", cex = 0.8)
        par(new = FALSE)
    }

    mtext(side = 3, plot.label, cex = 1.5)

    best.error <- sapply(tuned.list, function(x) x[which.min(x[,1]),1])
    best.exp <- sapply(tuned.list, function(x) x[which.max(x[,2]),2])

    error.mat <- cbind(min.error, best.error)
    colnames(error.mat) <- c("Maxnodes", "MSE")

    exp.mat <- cbind(max.exp, best.exp)
    colnames(exp.mat) <- c("Maxnodes", "Variance_Explained")

    result <- list(error.mat, exp.mat)
    invisible(result)
}