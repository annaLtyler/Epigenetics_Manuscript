#This function takes in a single response variable, and two different
#feature matrices. It then runs multiple random forests on each and 
#compares the mean squared error (MSE) for the two models across the
#range of parameters.
#to fix one of the factors, use a single number. To test one of the
#factors, use a range. Anything left as NULL will use the default 
#value.
#since nodesize and maxnodes constrain each other, we have chosen
#to include only maxnodes.
#to do a two dimensional parameter sweep, set ranges for both.
#The plots mark where the second model does better than the first.

compare.rf.models <- function(response, fm1, fm2, maxnodes = NULL, mtry = NULL, 
ntree = 500, verbose = FALSE, test.labels = c("fm1", "fm2"), filename = "results.RDS", 
plot.results = FALSE, plot.label = "", use.previous.results = TRUE, save.results = TRUE){

    if(is.null(mtry)){
        mtry1 <- floor(ncol(fm1)/3)
        mtry2 <- floor(ncol(fm2)/3)
    }else{
        mtry1 <- mtry
        mtry2 <- mtry
    }

    maxnodes.trials <- max(c(length(maxnodes), 1))
    mtry.trials <- max(c(length(mtry1), length(mtry2)))
    n.trials <- maxnodes.trials * mtry.trials


    if(!file.exists(filename) || !use.previous.results){
       
        if(verbose){cat("\nFitting", test.labels[1], "model", "\n")}

        mse.mat1 <- matrix(NA, nrow = maxnodes.trials, ncol = mtry.trials) 
        colnames(mse.mat1) <- paste0("mtry", mtry1)
        rownames(mse.mat1) <- paste0("maxnodes", maxnodes)

        idx = 1
        for(i in 1:mtry.trials){
            for(j in 1:maxnodes.trials){
                if(verbose){report.progress(idx, n.trials)}
                if(mtry1[i] <= ncol(fm1)){
                    full.rf <- randomForest(y = response, x = fm1, maxnodes = maxnodes[j], 
                    mtry = mtry1[i], ntree = ntree)
                    mse.mat1[j,i] <- median(full.rf$mse)
                }
                idx = idx + 1
            }
        }

        
        if(verbose){cat("\nFitting", test.labels[2], "model", "\n")}
        
        mse.mat2 <- matrix(NA, nrow = maxnodes.trials, ncol = mtry.trials) 
        colnames(mse.mat2) <- paste0("mtry", mtry2)
        rownames(mse.mat2) <- paste0("maxnodes", maxnodes)

        idx = 1
        for(i in 1:mtry.trials){
            for(j in 1:maxnodes.trials){
                if(verbose){report.progress(idx, n.trials)}
                if(mtry2[i] <= ncol(fm2)){
                    full.rf <- randomForest(y = response, x = fm2, maxnodes = maxnodes[j], 
                    mtry = mtry2[i], ntree = ntree)
                    mse.mat2[j,i] <- median(full.rf$mse)
                }
                idx = idx + 1
            }
        }


        if(verbose){cat("\n")}
    
        med.mse <- list(mse.mat1, mse.mat2)
        names(med.mse) <- test.labels
        if(save.results){
            saveRDS(med.mse, filename)
        }

    }else{
        med.mse <- readRDS(filename)
    }

    if(plot.results){

        if(maxnodes.trials == 1){
            plot.mat <- Reduce("rbind", med.mse)
            rownames(plot.mat) <- test.labels
        }
        if(mtry.trials == 1){
            plot.mat <- t(Reduce("cbind", med.mse))
            rownames(plot.mat) <- test.labels
        }
        
        if(mtry.trials == 1 || maxnodes.trials == 1){
            full.better <- which(plot.mat[1,] - plot.mat[2,] > 0)
            best <- which.max(plot.mat[1,] - plot.mat[2,])
            a <- barplot(plot.mat, beside = TRUE, legend.text = TRUE, 
            main = plot.label, las = 2)
            if(length(full.better) > 0){
                text(x = a[1,full.better], y = plot.mat[1,full.better], labels = "*", col = "blue")
                text(x = a[1,best], y = plot.mat[1,best], labels = "*", col = "red", cex = 1.5)
            }
        }else{
            pheatmap(med.mse[[1]], cluster_rows = FALSE, cluster_cols = FALSE, main = test.labels[1])
            pheatmap(med.mse[[2]], cluster_rows = FALSE, cluster_cols = FALSE, main = test.labels[2])
            pheatmap(med.mse[[1]] - med.mse[[2]], cluster_rows = FALSE, cluster_cols = FALSE, 
            main = paste(test.labels[1], "minus", test.labels[2]))
        }
    }

    invisible(med.mse)

}
