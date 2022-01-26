#This function plots multiple columns of a bed file.
#It is similar to plot_bed, but all columnes are 
#in a single file.

plot_combined_bed <- function(merged.bed.file, chr = 5, start.bp = 104459450, 
    end.bp = 104505819, bp.buffer = 1e6, col = "black", 
    value.start = 5, filename = "Peaks.pdf", label.names = NULL){

    bed.data <- read.table(merged.bed.file, stringsAsFactors = FALSE)

    chr.table <- bed.data[which(bed.data[,1] == chr),]
    plot.table <- chr.table[intersect(which(chr.table[,2] >= start.bp-bp.buffer), which(chr.table[,3] <= end.bp+bp.buffer)),]
    n.lines <- nrow(plot.table)

    if(n.lines == 0){cat("No data in this range\n");return(NULL)}

    just.vals <- plot.table[,value.start:ncol(plot.table)]
    if(is.null(label.names)){label.names = colnames(just.vals)}


    if(length(col) < ncol(just.vals)){col <- rep(col, ncol(just.vals))}
    min.x <- min(plot.table[,2])
    max.x <- max(plot.table[,3])
    min.y <- 0
    max.y <- sqrt(sqrt(max(just.vals)))


    pdf(filename)
    par(mar = c(1,6,0,0), xpd = NA)
    layout(matrix(1:ncol(just.vals), ncol = 1))
    for(i in 1:ncol(just.vals)){
        plot.new()
        plot.window(xlim = c(min.x, max.x), ylim = c(min.y, max.y))
        text(x = min.x, y = mean(c(min.y, max.y)), labels = label.names[i], 
            cex = 0.5, adj = 1)
        for(j in 1:nrow(plot.table)){
            draw.rectangle(min.x = plot.table[j,2], max.x = plot.table[j,3], 
            min.y = 0, max.y = sqrt(sqrt(just.vals[j,i])), fill = col[i], border = NA)
        }
    }
    segments(x0 = start.bp, x1 = end.bp, y0 = 0, lwd = 3)
    par(xpd = TRUE)
    dev.off()

    invisible(plot.table)



}