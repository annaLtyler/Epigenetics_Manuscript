#This function plots a methylation matrix with the 
#marks in their genomic position

plot.methyl.mat <- function(methyl.mat, plot.label = "", line.color = "black", 
xlim = NULL, ylab = "Methylation", bins = NULL){

    if(all(is.na(methyl.mat))){
        plot.text(paste("No methylation for", plot.label))
        return(NULL)
    }else{
        strain.cols <- matrix(colors.from.values(methyl.mat, use.pheatmap.colors = TRUE, 
        global.color.scale = TRUE, global.min = 0, global.max = 100), 
        ncol = ncol(methyl.mat), nrow = nrow(methyl.mat))

        methyl.pos <- as.numeric(colnames(methyl.mat))
        if(is.null(xlim)){xlim = c(min(methyl.pos), max(methyl.pos))}

        #trim the methylation matrix so it doesn't plot outside the plotting region
        use.idx <- intersect(which(methyl.pos >= xlim[1]), which(methyl.pos <= xlim[2]))
        methyl.mat <- methyl.mat[,use.idx]
        methyl.pos <- methyl.pos[use.idx]

        plot.new()
        plot.window(xlim = xlim, ylim = c(0.8, nrow(methyl.mat)+0.2))
        plot.dim <- par("usr")
        plot.width = plot.dim[2] - plot.dim[1]        

        par(xpd = TRUE)
        if(!is.null(bins)){
            num.islands <- max(bins)
            for(i in 1:num.islands){
                island.pos <- methyl.pos[which(bins == i)]
                draw.rectangle(min(island.pos), max(island.pos), 0, -0.5, fill = "gray")
            }
        }
        
        for(i in 1:nrow(methyl.mat)){
            points(x = methyl.pos, y = rep(i, length(methyl.pos)), 
            pch = "|", col = strain.cols[i,], cex = 1.5)
            text(plot.dim[1], i, rownames(methyl.mat)[i])
        }
        par(xpd = FALSE)
        mtext(ylab, side = 2)
        mtext(plot.label, side = 3)
    
        plot.dim <- par("usr")
        plot.width <- plot.dim[2] - plot.dim[1]
        plot.height <- plot.dim[4] - plot.dim[3]

        xmin <- plot.dim[2] + plot.width * 0.02
        xmax <- plot.dim[2] + plot.width * 0.06
        xmid <- mean(c(xmin, xmax))
        ymin <- plot.dim[3]
        ymax <- plot.dim[4]
        state.cols <- colors.from.values(c(0, 50, 100), use.pheatmap.colors = TRUE)
        state.labels <- c(0, 50, 100)
        yseg <- segment.region(ymin, ymax, 4, alignment = "ends")
        par(xpd = TRUE)
        for(i in 1:(length(yseg)-1)){
            draw.rectangle(xmin, xmax, yseg[i], yseg[i+1], fill = state.cols[i], 
            border.col = line.color)
            text(x = xmid, y = mean(c(yseg[i], yseg[i+1])), labels = state.labels[i])
        }

        par(xpd = FALSE)

    }


}