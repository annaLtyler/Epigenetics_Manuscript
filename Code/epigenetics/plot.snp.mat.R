#This function is similar to plot.chrom.mat, but doesn't 
#fill in regions between SNPs

plot.snp.mat <- function(snp.mat, num.states = 4, xlim = NULL, 
state.labels = 1:num.states, ylab = "SNP", line.color = "gray", 
state.cols = NULL, empty.cell.color = "lightgray"){

    #pheatmap(state.mat, cluster_rows = FALSE, cluster_cols = FALSE)
    if(is.null(state.cols)){
        state.cols <- colors.from.values(1:num.states, use.pheatmap.colors = TRUE, 
        global.color.scale = TRUE, global.min = 1, global.max = num.states)
    }

    snp.pos <- as.numeric(colnames(snp.mat))
    if(is.null(xlim)){xlim <- c(min(snp.pos), max(snp.pos))}
  
    #quartz()
    plot.new()
    plot.window(xlim = xlim, ylim = c(0.8, nrow(snp.mat)+0.2))
    plot.dim <- par("usr")
    plot.width <- plot.dim[2] - plot.dim[1]
    plot.height <- plot.dim[4] - plot.dim[3]
    #draw.rectangle(plot.dim[1], plot.dim[2], plot.dim[3], plot.dim[4])

    par(xpd = TRUE)
    for(i in 1:nrow(snp.mat)){
        ypos <- nrow(snp.mat) - i + 1
        text(plot.dim[1], ypos, rownames(snp.mat)[i])
        seg.col <- rep(empty.cell.color, ncol(snp.mat))
        for(sn in 1:num.states){
            seg.col[which(snp.mat[i,] == sn)] <- state.cols[sn]
        }
        segments(x0 = snp.pos, x1 = snp.pos, y0 = ypos-0.5, y1 = ypos+0.5, lwd = 2,
            col = seg.col)
    
    }
    par(xpd = FALSE)
    mtext(ylab, side = 2)

    xmin <- plot.dim[2] + plot.width * 0.02
    xmax <- plot.dim[2] + plot.width * 0.06
    xmid <- mean(c(xmin, xmax))
    ymin <- plot.dim[3]
    ymax <- plot.dim[4]
    yseg <- segment.region(ymax, ymin, num.states+1, alignment = "ends")
    par(xpd = TRUE)
    for(i in 1:(length(yseg)-1)){
        draw.rectangle(xmin, xmax, yseg[i], yseg[i+1], fill = state.cols[i], 
        border.col = line.color)
        text(x = xmid, y = mean(c(yseg[i], yseg[i+1])), labels = state.labels[i])
    }

    par(xpd = FALSE)


    
}
