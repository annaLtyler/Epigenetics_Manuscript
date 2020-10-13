#This is a test function to see if we can find really
#good windows between methylation and expression.

plot.best.methyl.cor <- function(gene.name, methyl.mat, ordered.expr, gene.start,
gene.end, strand, min.window = 2, max.window = 200, plot.results = FALSE){

    if(max.window > ncol(methyl.mat)){max.window <- round(ncol(methyl.mat)/2)}
    
    get.coef <- function(test, type = c("beta", "sd", "t")){
        coef.mat <- summary(test)$coefficients
        if(nrow(coef.mat) == 2){
            if(type == "beta"){
                result = coef.mat[2,"Estimate"]
            }
            if(type == "sd"){
                result <- coef.mat[2,"Std. Error"]
            }
            if(type == "t"){
                result <- coef.mat[2,"t value"]
            }
        }else{
            result <- NA
        }
        return(result)
    }


    test.window <- function(ordered.expr, methyl.x){
        result <- try(lm(as.vector(ordered.expr)~as.vector(methyl.x)))
        if(class(result) == "try-class error"){
            return(NA)
            }else{
                return(result)
            }
    }

    window.v <- min.window:max.window
    all.window.effect <- all.window.x <- vector(mode = "list", length = length(window.v))
    names(all.window.effect) <- names(all.window.x) <- min.window:max.window
    idx <- 1
    for(i in window.v){
        window.idx <- sliding.window.el(1:ncol(methyl.mat), i, 1)
        mean.x <- sapply(window.idx, function(x) mean(as.numeric(colnames(methyl.mat)[x])))
        all.window.x[[idx]] <- mean.x
        mean.methyl <- lapply(window.idx, function(x) apply(methyl.mat[,x], 1, function(y) mean(y, na.rm = TRUE)))
        #all.window.effect[[idx]] <- round(suppressWarnings(sapply(mean.methyl, function(x) cor(x, ordered.expr))), 2)
        dummy.methyl <- lapply(mean.methyl, function(x) matrix(rep(x, nrow(ordered.expr)), nrow = nrow(ordered.expr), byrow = TRUE ))
        window.test <- lapply(dummy.methyl, function(x) test.window(ordered.expr, x))
        window.beta <- sapply(window.test, function(x) get.coef(x, "beta"))
        window.sd <- sapply(window.test, function(x) get.coef(x, "sd"))
        window.t <- sapply(window.test, function(x) get.coef(x, "t"))
        results.table <- rbind(window.beta, window.sd, window.t)
        colnames(results.table) <- mean.x
        all.window.effect[[idx]] <- results.table
        idx <- idx + 1
    }


    if(plot.results){
        #look at the position of the correlations
        cols <- colors.from.values(window.v, use.pheatmap.colors = TRUE)
        #cols <- colors.from.values(window.v, col.scale = "purple", light.dark = "d")

        #quartz(width = 10, height = 5)
        layout(matrix(c(1,2), nrow = 1), widths = c(1,0.2))
        plot.new()
        xlim <- c(min(as.numeric(colnames(methyl.mat))), max(as.numeric(colnames(methyl.mat))))
        ymin <- min(c(min(unlist(lapply(all.window.effect, function(x) x[1,])), na.rm = TRUE), -1))
        ymax <- max(c(max(unlist(lapply(all.window.effect, function(x) x[1,])), na.rm = TRUE), 1))
        plot.window(xlim = xlim, ylim = c(ymin, ymax))
        for(i in 1:length(all.window.effect)){
            points(all.window.x[[i]], all.window.effect[[i]][1,], type = "l", col = cols[i])
        }
        axis(1); axis(2)
        abline(h = 0)
        mtext(paste("Effect of Methylation on", gene.name, "Expression"))
        par(xpd = TRUE)
        if(strand == 1){
            arrows(x0 = gene.start, x1 = gene.end, y0 = ymax, y1 = ymax, lwd = 2)
        }else{
            arrows(x0 = gene.end, x1 = gene.start, y0 = ymax, y1 = ymax, lwd = 2)
        }
        par(xpd = FALSE)
        mtext("Effect Size (beta)", side = 2, line = 2)
        mtext("Mean Genomic Position of Window", side = 1, line = 2)

        imageWithTextColorbar(matrix(window.v, ncol = 1), use.pheatmap.colors = TRUE, cex = 1)
        mtext("Window Size", side = 3, line = -2)

        layout(matrix(c(1), nrow = 1))
        par(mar = c(5.1, 4.1, 4.1, 2.1))
    }
    
    invisible(all.window.effect)

}
