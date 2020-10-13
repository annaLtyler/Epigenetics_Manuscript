#This function plots a methylation matrix with the 
#marks in their genomic position

plot.methyl.mat <- function(methyl.mat, plot.label = "", xlim = NULL, bins = NULL){

    if(all(is.na(methyl.mat))){
        plot.text(paste("No methylation for", plot.label))
        return(NULL)
    }else{
        strain.cols <- matrix(colors.from.values(methyl.mat, use.pheatmap.colors = TRUE, 
        global.color.scale = TRUE, global.min = 0, global.max = 100), 
        ncol = ncol(methyl.mat), nrow = nrow(methyl.mat))

        methyl.pos <- as.numeric(colnames(methyl.mat))
        if(is.null(xlim)){xlim = c(min(methyl.pos), max(methyl.pos))}

        plot.new()
        plot.window(xlim = xlim, ylim = c(0,nrow(methyl.mat)+1))
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
            text(plot.dim[2], i, rownames(methyl.mat)[i])
        }
        par(xpd = FALSE)
        mtext("Methylation", side = 2)
        mtext(plot.label, side = 3)
    }
}