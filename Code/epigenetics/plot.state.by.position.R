#This function plots the results of one state.by.position()

plot.state.by.position <- function(state.pos.results, error.type = c("sd", "se", "none"),
plot.label = "", col = "gray", x.axis = TRUE, xlim = NULL, ylim = NULL, add = FALSE, lwd = 2){

    error.type = error.type[1]

    pos <- state.pos.results$position
    state.avg <- state.pos.results$state.avg

    if(error.type == "sd"){
        state.error <- state.pos.results$state.sd
    }
    if(error.type == "se"){
        state.error <- state.pos.results$state.se
    }
    if(error.type == "none"){
        state.error <- rep(0, length(state.avg))
    }

    if(is.null(ylim)){
        ymax <- max((state.avg + state.error), na.rm = TRUE)
        ymin <- min((state.avg - state.error), na.rm = TRUE)
        ylim <- c(ymin, ymax)
    }

    if(is.null(xlim)){
        xmin <- min(pos)
        xmax <- max(pos)
        xlim <- c(xmin, xmax)
    }

    if(!add){
        plot(pos, state.avg, main = plot.label, type = "l", xlab = "Relative Gene Position", 
        ylab = paste("Average State Presence"), pch = 16, cex = 0.7, col = col,
        ylim = ylim, axes = FALSE, xlim = xlim, lwd = lwd)
    }else{
        points(pos, state.avg, pch = 16, cex = 0.7, col = col, ylim = c(ymin, ymax), type = "l",
        lwd = lwd)
    }
    axis(2)
    if(x.axis){axis(1)}
    error.col <- col2rgb(col)
    plot.poly.xy(poly.top.x = pos, poly.top.y = state.avg + state.error, 
    poly.bottom.x = pos, poly.bottom.y = state.avg - state.error, 
    col = rgb(error.col[1]/256, error.col[2]/256, error.col[3]/256, alpha = 0.5), border = NA)
    abline(v = c(0,1), col = "darkgray")

}