#This function compares the importance measures of a randomForest object
#to see if there are any major disagreements.

compare.importance <- function(rf, transform = TRUE){

    import <- importance(rf)

    if(transform){
        import <- apply(import, 2, rankZ)
        xlab <- paste("Transformed", colnames(import)[1])
        ylab <- paste("Transformed", colnames(import)[2])
    }else{
        xlab <- colnames(import)[1]
        ylab <- colnames(import)[2]
    }

    min.x <- min(import[,1])
    max.x <- max(import[,1])
    plot.width <- max.x - min.x

    min.y <- min(import[,2])
    max.y <- max(import[,2])
    plot.height <- max.y - min.y

    plot.with.model(import[,1], import[,2], xlab = xlab, ylab = ylab, 
    ylim = c(min.y, max.y), xlim = c(min.x - (plot.width*0.2), max.x), 
    cex = 0.7, pch = 16)
    text(import[,1], import[,2], rownames(import), pos = 2, offset = 0.2)
}