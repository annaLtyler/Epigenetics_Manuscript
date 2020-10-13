tune.maxnodes <- function(y, x, maxnodes = 1:100, mtry = NULL, ntree = 500, 
filename = "Tuned.RF.RDS", plot.results = TRUE, plot.label = "", 
use.previous.results = TRUE, save.results = TRUE, verbose = FALSE, 
keep.forest = FALSE){
    
    if(is.null(mtry)){
        mtry <- floor(ncol(x)/3)
    }else{
        mtry <- mtry
    }

    if(!file.exists(filename) || !use.previous.results){

        n.trials <- length(maxnodes)
        results.mat <- matrix(NA, nrow = n.trials, ncol = 2)
        colnames(results.mat) <- c("Med_MSE", "Var_Exp")
        for(i in 1:n.trials){
            if(verbose){report.progress(i, n.trials)}
            full.rf <- randomForest(y = y, x = x, maxnodes = maxnodes[i], 
            mtry = mtry, ntree = ntree, keep.forest = keep.forest)
            results.mat[i,1] <- tail(full.rf$mse, 1)
            results.mat[i,2] <- rf.var(full.rf)
        }
    if(save.results){
        saveRDS(results.mat, filename)
    }
    }else{
        results.mat <- readRDS(filename)
    }    

    if(plot.results){
        
        par(mar = c(5,5,5,5))
        plot.new()
        plot.window(xlim = c(min(maxnodes), max(maxnodes)), 
        ylim = c(min(results.mat[,1]), max(results.mat[,1])))
        line.cols <- c("#d8b365", "#5ab4ac")
        plot(maxnodes, results.mat[,1], ylab = "MSE", type = "l", 
        col = line.cols[1], lwd = 2)
        par(new = TRUE)
        plot(maxnodes, results.mat[,2], xlab = NA, ylab = NA, type = "l", 
        col = line.cols[2], lwd = 2, axes = FALSE)
        axis(side = 4)
        mtext(side = 4, line = 2.5, "Variance Explained")
        
        par(xpd = TRUE)
        plot.dim <- par("usr")
        plot.height <- plot.dim[4] - plot.dim[3]
        plot.width <- plot.dim[2] - plot.dim[1]
        legend(plot.dim[1], (plot.dim[4]+(plot.height*0.08)), 
        col = line.cols, lty = 1, legend = c("MSE", "Variance Explained"), 
        horiz = TRUE, lwd = 2)
        mtext(side = 3, line = 2.5, plot.label)

        par(xpd = FALSE)
        max.pt <- which.max(results.mat[,2] - results.mat[,1])
        abline(v = max.pt, lwd = 2, col = "darkgray", lty = 2)
        par(new = FALSE)
    }

    invisible(results.mat)

}