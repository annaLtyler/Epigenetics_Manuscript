#This function plots peaks in bedfiles at a specific
#location


plot_bed <- function(bed.files, chr = 5, start.bp = 104459450, 
    end.bp = 104505819, bp.buffer = 1e6, col = "black", 
    value.col = 4, filename = "Peaks.pdf"){

    if(length(bed.files) > length(col)){rep(col, length(bed.files))}

    bed.data <- lapply(bed.files, function(x) read.table(x, stringsAsFactors = FALSE))

    chr.tables <- lapply(bed.data, function(x) x[which(x[,1] == chr),])
    plot.tables <- lapply(chr.tables, function(x) x[intersect(which(x[,2] >= start.bp-bp.buffer), which(x[,3] <= end.bp+bp.buffer)),])
    n.lines <- sapply(plot.tables, nrow)

    if(all(n.lines == 0)){cat("No data in this range\n");return(NULL)}

    min.x <- min(sapply(plot.tables, function(x) if(nrow(x) > 0){min(x[,2])}else{NA}), na.rm = TRUE)
    max.x <- max(sapply(plot.tables, function(x) if(nrow(x) > 0){max(x[,3])}else{NA}), na.rm = TRUE)
    min.y <- 0
    max.y <- max(sapply(plot.tables, function(x) if(nrow(x) > 0){max(x[,value.col])}else{NA}), na.rm = TRUE)


    pdf(filename)
    par(mar = c(1,6,0,0), xpd = NA)
    layout(matrix(1:length(plot.tables), ncol = 1))
    for(i in 1:length(plot.tables)){
        plot.new()
        plot.window(xlim = c(min.x, max.x), ylim = c(min.y, max.y))
        if(nrow(plot.tables[[i]]) > 0){
            text(x = min.x, y = mean(c(min.y, max.y)), labels = basename(bed.files)[i], 
                cex = 0.5, adj = 1)
            for(j in 1:nrow(plot.tables[[i]])){
                draw.rectangle(min.x = plot.tables[[i]][j,2], max.x = plot.tables[[i]][j,3], 
                min.y = 0, max.y = plot.tables[[i]][j,value.col], fill = col[i], border = NA)
            }
        }
    }
    segments(x0 = start.bp, x1 = end.bp, y0 = 0, lwd = 3)
    par(xpd = TRUE)
    dev.off()

    invisible(plot.tables)

}

