#This function is similar to plot.state.mat, but plots a
#chromatin matrix on its own using the bp positions and 
#state colors
#island.bins is for drawing outlines of methyl bins

plot.chrom.mat <- function(state.mat, num.states = 8, xlim = NULL, 
state.labels = 1:num.states, ylab = "Chromatin State", island.bins = NULL, 
line.color = "gray", state.cols = NULL, empty.cell.color = "lightgray"){

    #pheatmap(state.mat, cluster_rows = FALSE, cluster_cols = FALSE)
    if(is.null(state.cols)){
        state.cols <- colors.from.values(1:num.states, use.pheatmap.colors = TRUE, 
        global.color.scale = TRUE, global.min = 1, global.max = num.states)
    }

    methyl.pos <- as.numeric(colnames(state.mat))
    if(is.null(xlim)){xlim <- c(min(methyl.pos, na.rm = TRUE), max(methyl.pos, na.rm = TRUE))}
  
    #quartz()
    plot.new()
    plot.window(xlim = xlim, ylim = c(0.8, nrow(state.mat)+0.2))
    plot.dim <- par("usr")
    plot.width <- plot.dim[2] - plot.dim[1]
    plot.height <- plot.dim[4] - plot.dim[3]
    #draw.rectangle(plot.dim[1], plot.dim[2], plot.dim[3], plot.dim[4])

    par(xpd = TRUE)
    for(i in 1:nrow(state.mat)){
        ypos <- nrow(state.mat) - i + 1
        text(plot.dim[1], ypos, rownames(state.mat)[i])
        consec.states <- consec.pairs(state.mat[i,])
        transition.pts <- which(apply(consec.states, 1, function(x) x[1] != x[2]))

        #get the first state in the row
        first.state <- consec.states[1,1]
        first.state.start <- methyl.pos[1]
        state.col <- state.cols[first.state]
        if(is.na(state.col)){state.col <- empty.cell.color}

        #if there are no changes in state, use the first state for the whole row
        if(length(transition.pts) == 0){    
            draw.rectangle(min(methyl.pos), max(methyl.pos), ypos-0.5, 
            ypos+0.5, fill = state.col, border = NA)
        }else{ #otherwise, draw
            for(tp in transition.pts){
                next.state <- consec.states[tp,2]
                next.state.start <- methyl.pos[tp+1]
                next.state.col <- state.cols[next.state]
                draw.rectangle(first.state.start, next.state.start, ypos-0.5, 
                ypos+0.5, fill = state.col, border = NA)
                first.state.start <- next.state.start
                state.col <- next.state.col
            }
            #add the last section
            next.state.start <- tail(methyl.pos, 1)
            draw.rectangle(first.state.start, next.state.start, ypos-0.5, 
            ypos+0.5, fill = state.col, border = NA)
        }
    }
    par(xpd = FALSE)
    mtext(ylab, side = 2)

    if(!is.null(island.bins)){
        u_islands <- unique(island.bins[,1])
        for(i in 1:length(u_islands)){
            island.locale <- which(island.bins[,1] == u_islands[i])
            island.start <- min(island.bins[island.locale,2])
            island.end <- max(island.bins[island.locale,2])
            draw.rectangle(island.start, island.end, 0.5, nrow(state.mat)+0.5, border = line.color)
        }
    }


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