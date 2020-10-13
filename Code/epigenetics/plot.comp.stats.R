#This function plots a series of results obtained from compare2tuned.
#See Explore_RF.Rmd for an example as to how to build a list with
#compare2tuned.R


plot.comp.stats <- function(comp.stats, plot.type = c("boxplot", "plot", "by.trial"), 
plot.title = "", comp.labels = NULL, trial.labels = NULL){

    par(mar = c(4,4,4,1))
    plot.type <- plot.type[1]

    if(is.null(comp.labels)){
        comp.labels <- rownames(comp.stats[[1]][[1]])
    }


  sim.error <- t(sapply(comp.stats, function(x) x[[1]][1,]))
  gen.error <- t(sapply(comp.stats, function(x) x[[1]][2,]))

  sim.exp <- t(sapply(comp.stats, function(x) x[[2]][1,]))
  gen.exp <- t(sapply(comp.stats, function(x) x[[2]][2,]))

    

    if(plot.type == "boxplot"){
        par(mfrow = c(2,2))
        boxplot(list(sim.error[,2], gen.error[,2]), main = "Lowest MSE", 
        names = comp.labels)

        boxplot(list(sim.exp[,2], gen.exp[,2]), main = "Highest Variance Explained", 
        names = comp.labels)

        boxplot(list(sim.error[,1], gen.error[,1]), main = "Nodes Used for Lowest MSE", 
        names = comp.labels)

        boxplot(list(sim.exp[,1], gen.exp[,1]), 
        main = "Nodes Used for Highest\nVariance Explained", names = comp.labels)
    }

    if(plot.type == "plot"){
        par(mfrow = c(2,2))
        plot(sim.error[,2], gen.error[,2], xlab = paste("MSE for", comp.labels[1]), 
        ylab = paste("MSE for", comp.labels[2]), main = "MSE")
        abline(0,1)

        plot(sim.exp[,2], gen.exp[,2] , 
        xlab = paste("Variance Explained for", comp.labels[1]), 
        ylab = paste("Variance Explained for", comp.labels[2]), 
        main = "Variance Explained")
        abline(0,1)

        plot(sim.error[,1], gen.error[,1], xlab = paste("Maxnodes for", comp.labels[1]), 
        ylab = paste("Maxnodes for", comp.labels[2]), main = "Maxnodes for MSE")
        abline(0,1)

        plot(sim.exp[,1], gen.exp[,1], xlab = paste("Maxnodes for", comp.labels[1]), 
        ylab = paste("Maxnodes for", comp.labels[2]), 
        main = "Maxnodes for Variance Explained")
        abline(0,1)

    }

    if(plot.type == "by.trial"){
        cols <- c("#d8b365", "#5ab4ac")
        
        layout.mat <- matrix(c(1,2,5, 3,4,0), nrow = 2, byrow = TRUE)
        layout(layout.mat, widths = c(1,1,0.3))
        if(is.null(trial.labels)){trial.labels <- 1:nrow(gen.error)}
        
        barplot(rbind(sim.error[,2], gen.error[,2]), names = trial.labels, 
        main = "Lowest MSE", beside = TRUE, col = cols)
        
        barplot(rbind(sim.exp[,2], gen.exp[,2]), names = trial.labels, 
        main = "Highest Variance Explained", beside = TRUE, col = cols)
        
        barplot(rbind(sim.error[,1], gen.error[,1]), names = trial.labels, 
        main = "Nodes Used for Lowest MSE", beside = TRUE, col = cols)
        
        barplot(rbind(sim.exp[,1], gen.exp[,1]), names = trial.labels, 
        main = "Nodes Used for Highest Variance Explained", beside = TRUE, 
        col = cols)

        par(mar = c(0,0,0,0))
        plot.new()
        plot.window(xlim = c(0,1), ylim = c(0,1))
        legend(x = 0, y = 0.9, legend = comp.labels, fill = cols)
    }
    mtext(plot.title, side = 3, outer = TRUE, line = -1.5)
}