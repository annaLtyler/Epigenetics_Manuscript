#This barplotting function plots numbers on
#top of each bar.


barplot_with_num <- function(values, xlab = "", ylab = "", main = "", col = "gray",
    cex = 1, text.gap = 0.1, text.shift = 0.05, horiz = FALSE, names = NULL, las = 1,
    srt = 0, adj = 0.5){

    plot.range <- max(values, na.rm = TRUE) - min(values, na.rm = TRUE)
    coords <- barplot(values, plot = FALSE)

    if(!horiz){
        ylim <- c(0, max(values, na.rm = TRUE)+(plot.range*text.gap))
        xlim <- c(0, max(coords, na.rm = TRUE))
    }else{
        ylim <- c(0, max(coords, na.rm = TRUE))
        xlim <- c(0, max(values, na.rm = TRUE)+(plot.range*text.gap))
    }

    a <- barplot(values, xlab = xlab, ylab = ylab, main = main, col = col, 
        ylim = ylim, xlim = xlim, horiz = horiz, names = names, las = las)
    text(a[,1], y = values+(plot.range*text.shift), labels = values, cex = cex, 
        srt = srt, adj = adj)

}