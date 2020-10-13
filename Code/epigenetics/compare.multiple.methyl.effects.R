#This function takes in a list of results from 
#methyl.expression.cor and looks for common areas
#of associations between methylation and expression
#across multiple genes. 
#we break each gene into multiple sections: upstream, 
#downstream, gene body, and first exon.

compare.multiple.methyl.effects <- function(methyl.results.list, window.size = 20, max.bp = 1000){

    #====================================================================
    # internal functions
    #====================================================================
    get.flanking <- function(methyl.result, side = c("up", "down")){
        gene.start <- methyl.result$gene.info$gene_start
        gene.end <- methyl.result$gene.info$gene_end
        gene.strand <- methyl.result$gene.info$gene_strand

        #This code is for getting upstream sequence
        #if we want downstream, flip the strand
        if(side == "down"){
            gene.strand <- gene.strand * -1
        }

        window.locale <- which(names(methyl.result$window.cor) == window.size)
        window.result <- methyl.result$window.cor[[window.locale]]

        if(gene.strand == 1){
            upstream.locale <- which(as.numeric(colnames(window.result)) < gene.start)
        }else{
            upstream.locale <- which(as.numeric(colnames(window.result)) > gene.end)
        }
        upstream.windows <- window.result[1,upstream.locale]
        pos <- as.numeric(colnames(window.result)[upstream.locale])
        rel.pos <- pos - min(pos)
        names(upstream.windows) <- rel.pos
        #plot(upstream.windows, type = "l")
        return(upstream.windows)
    }

    get.first.exon <- function(methyl.result){
        gene.start <- methyl.result$gene.info$gene_start
        gene.end <- methyl.result$gene.info$gene_end
        gene.strand <- methyl.result$gene.info$gene_strand
        exon.start <- methyl.result$gene.info$exon_start
        exon.end <- methyl.result$gene.info$exon_end

        if(gene.strand == 1){
            near.idx <- get.nearest.pt(exon.start, gene.start)
            #first.start <- min(exon.start)
            #first.idx <- which.min(exon.start)
            #first.end <- exon.end[first.idx]
        }else{
            near.idx <- get.nearest.pt(exon.start, gene.end)
            #first.start <- max(exon.start)
            #first.idx <- which.max(exon.start)
            #first.end <- exon.end[first.idx]
        }
        first.start <- exon.start[near.idx]
        first.end <- exon.end[near.idx]

        window.locale <- which(names(methyl.result$window.cor) == window.size)
        window.result <- methyl.result$window.cor[[window.locale]]

        methyl.pos <- as.numeric(colnames(window.result))
        exon.locale <- intersect(which(methyl.pos >= first.start), which(methyl.pos <= first.end))
        exon.window <- window.result[1,exon.locale]
        pos <- as.numeric(colnames(window.result)[exon.locale])
        rel.pos <- pos - min(pos)
        names(exon.window) <- rel.pos
        return(exon.window)
    }

    get.gene.body <- function(methyl.result){
        gene.start <- methyl.result$gene.info$gene_start
        gene.end <- methyl.result$gene.info$gene_end
        gene.strand <- methyl.result$gene.info$gene_strand

        window.locale <- which(names(methyl.result$window.cor) == window.size)
        window.result <- methyl.result$window.cor[[window.locale]]

        methyl.pos <- as.numeric(colnames(window.result))
        gene.locale <- intersect(which(methyl.pos >= gene.start), which(methyl.pos <= gene.end))
        gene.window <- window.result[1,gene.locale]
        pos <- as.numeric(colnames(window.result)[gene.locale])

        if(gene.strand == -1){
            rel.pos <- -pos - min(-pos)
        }else{
            rel.pos <- pos - min(pos)
        }
        names(gene.window) <- rel.pos
        return(gene.window)

    }


    plot.aligned <- function(window.list){

        all.pos <- unlist(lapply(window.list, function(x) as.numeric(names(x))))
        max.windows <- max(all.pos)
        min.windows <- 0

        all.vals <- unlist(window.list)
        min.val <- min(all.vals, na.rm = TRUE)
        max.val <- max(all.vals, na.rm = TRUE)

        layout(matrix(c(1,2), nrow = 1), widths = c(1,0.3))
        plot.new()
        plot.window(xlim = c(min.windows, max.windows), ylim = c(0, length(window.list)))
        for(i in 1:length(window.list)){
            if(length(window.list[[i]]) > 1){
                pt.cols <- colors.from.values(window.list[[i]], use.pheatmap.colors = TRUE,
                global.color.scale = TRUE, global.min = min.val, global.max = max.val)
                points(x = as.numeric(names(window.list[[i]])), 
                y = rep(i, length(window.list[[i]])), col = pt.cols, pch = "|", cex = 1.5)
            }
        }
        imageWithTextColorbar(matrix(seq(min.val, max.val, 0.1), ncol = 1), 
        use.pheatmap.colors = TRUE, cex = 1)
    }


    #====================================================================
    # end internal functions
    #====================================================================


    not.null <- which(!sapply(methyl.results.list, is.null))
    test.down <- lapply(methyl.results.list[not.null], function(x) get.flanking(x, "down"))
    test.up <- lapply(methyl.results.list[not.null], function(x) get.flanking(x, "up"))
    test.exon <- lapply(methyl.results.list[not.null], get.first.exon)
    test.gene <- lapply(methyl.results.list[not.null], get.gene.body)

    plot.aligned(test.down)
    plot.aligned(test.up)
    plot.aligned(test.exon)
    plot.aligned(test.gene)
    
    boxplot(test.up);abline(h = 0)
    boxplot(test.down);abline(h = 0)
    boxplot(test.exon);abline(h = 0)
    boxplot(test.gene);abline(h = 0)



}