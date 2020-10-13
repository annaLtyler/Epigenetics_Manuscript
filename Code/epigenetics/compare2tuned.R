#This function takes in two objects returned from tune.maxnodes
#and plots them on the same axes to compare.

compare2tuned <- function(tuned1, tuned2, plot.labels = c("Full", "Genetic"), 
plot.results = TRUE, legend.y = 0.17){

    xmin <- 0; xmax <- nrow(tuned1)
    ymin1 <- min(c(tuned1[,1], tuned2[,1]))
    ymin2 <- min(c(tuned1[,2], tuned2[,2]))
    ymax1 <- max(c(tuned1[,1], tuned2[,1]))
    ymax2 <- max(c(tuned1[,2], tuned2[,2]))

    min.error1 <- which.min(tuned1[,1])
    min.error2 <- which.min(tuned2[,1])
    max.exp1 <- which.max(tuned1[,2])
    max.exp2 <- which.max(tuned2[,2])

    if(plot.results){
        line.cols <- c("#a6611a", "#d8b365", "#018571", "#5ab4ac")

        par(mar = c(5,5,5,5))
        plot.new()
        plot.window(xlim = c(xmin, xmax), ylim = c(ymin1, ymax1))
        
        plot(1:nrow(tuned1), tuned1[,1], ylab = "MSE", type = "l", 
        col = line.cols[1], lwd = 2, xlab = "Number of Terminal Nodes",
        ylim = c(ymin1, ymax1))
        points(1:nrow(tuned2), tuned2[,1], type = "l", col = line.cols[2], lwd = 2)
        points(x = c(min.error1, min.error2), y = c(tuned1[min.error1,1], tuned2[min.error2,1]),
        pch = "*", col = line.cols[1:2], cex = 2.5)
    
        par(new = TRUE)
        plot(1:nrow(tuned1), tuned1[,2], xlab = NA, ylab = NA, type = "l", 
        col = line.cols[3], lwd = 2, axes = FALSE, ylim = c(ymin2, ymax2))
        points(1:nrow(tuned2), tuned2[,2], type = "l", col = line.cols[4], lwd = 2)
        axis(side = 4)
        mtext(side = 4, line = 2.5, "Variance Explained")
        points(x = c(max.exp1, max.exp2), y = c(tuned1[max.exp1,2], tuned2[max.exp2,2]),
        pch = "*", col = line.cols[3:4], cex = 2.5)
    

        par(xpd = TRUE)
        plot.dim <- par("usr")
        plot.height <- plot.dim[4] - plot.dim[3]
        plot.width <- plot.dim[2] - plot.dim[1]

        legend(x = plot.dim[1], y = (plot.dim[4]+(plot.height*legend.y)), 
        col = line.cols[1:2], lty = 1, legend = c(plot.labels[1], plot.labels[2]), 
        horiz = FALSE, lwd = 2,title = "MSE")

        legend(x = (plot.dim[1]+(plot.width/2)), y = (plot.dim[4]+(plot.height*legend.y)), 
        col = line.cols[3:4], lty = 1, legend = c(plot.labels[1], plot.labels[2]), 
        horiz = FALSE, lwd = 2,title = "Variance Explained")
        par(new = FALSE)
    }

    min.error1 <- which.min(tuned1[,1])
    min.error2 <- which.min(tuned2[,1])
    max.exp1 <- which.max(tuned1[,2])
    max.exp2 <- which.max(tuned2[,2])

    error.mat <- matrix(c(min.error1, tuned1[min.error1,1], min.error2, tuned2[min.error2,1]),
    nrow = 2, byrow = TRUE)
    colnames(error.mat) <- c("Maxnodes", "MSE")
    rownames(error.mat) <- plot.labels

    exp.mat <- matrix(c(max.exp1, tuned1[max.exp1,2], max.exp2, tuned2[max.exp2,2]),
    nrow = 2, byrow = TRUE)
    colnames(exp.mat) <- c("Maxnodes", "Variance_Explained")
    rownames(exp.mat) <- plot.labels

    result <- list(error.mat, exp.mat)
    invisible(result)
}