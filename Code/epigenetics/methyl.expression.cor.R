#This function finds the correlation between 
#gene expression and DNA methylation for a
#named gene
#gene.frac.gap defines what kind of gap between
#methyl marks is considered big enough to create
#a new island. The default is 1% of the gene length

methyl.expression.cor <- function(gene.name, upstream.buffer = 5000, 
downstream.buffer = 5000, rrbs.data, gene.info, rna.seq, strains, gene.frac.gap = 0.01,
min.representation = 3, treatment = "C", col.table, id.table, sliding.window = FALSE, 
plot.results = FALSE, save.images = FALSE, path = ".", plot.label = ""){

    gene.locale <- which(gene.info[,"external_gene_name"] == gene.name)
    if(length(gene.locale) == 0){stop("can't find ", gene.name)}

    gene.chr <- unique(gene.info[gene.locale,"chromosome_name"])
    gene.start <- unique(gene.info[gene.locale,"start_position"])
    gene.end <- unique(gene.info[gene.locale,"end_position"])
    gene.strand <- unique(gene.info[gene.locale,"strand"])
    exon.start <- unique(gene.info[gene.locale,"exon_chrom_start"])
    exon.end <- unique(gene.info[gene.locale,"exon_chrom_end"])

    gene.info.table <- list("gene.name" = gene.name, "Chr" = gene.chr, 
    "gene_start" = gene.start, "gene_end" = gene.end, "gene_strand" = gene.strand, 
    "exon_start" = exon.start, "exon_end" = exon.end)

    gene.expr <- get.gene.expr(gene.info, rna.seq, gene.name, 
    strain.name = strains, treatment = treatment, col.table = col.table, 
    average.replicates = FALSE)
    colnames(gene.expr) <- substr(colnames(gene.expr), 1, 2)
    mean.expr <- colMeans(gene.expr)
    expr.order <- order(mean.expr)

    if(plot.results){
        if(save.images){
            jpeg(file.path(path, paste0(plot.label, ".Expression.", gene.name, ".jpg")),
            width = 6, height = 5, units = "in", res = 300)
        }else{
                quartz(width = 6, height = 5)
        }
        boxplot(gene.expr[,expr.order], 
        main = paste("Expression for", gene.name, "across strains"))
        if(save.images){
            dev.off()
        }
    }

    methyl.table <- get.methyl.pos(all.methyl.data = rrbs.data, 
    strain.name = strains, treatment.name = treatment, 
    id.table, chr = gene.chr, start.pos = gene.start - upstream.buffer, 
    stop.pos = gene.end + downstream.buffer, average.replicates = TRUE)
    if(class(methyl.table) == "character"){
        final.results <- list("gene.info" = gene.info.table, "expr" = gene.expr, 
        "methylation.matrix" = methyl.table)
        return(final.results)
    }
    colnames(methyl.table) <- gsub(paste0("_", treatment), "", colnames(methyl.table))

    if(length(strains) > 1){
        methyl.order <- order.strains(colnames(gene.expr[,expr.order,drop=FALSE]), 
        colnames(methyl.table), col.table)
        methyl.order <- methyl.order[which(!is.na(methyl.order))]
    }else{
        methyl.order <- which(colnames(methyl.table) == strains)
    }
    

    methyl.mat <- t(methyl.table[,methyl.order]) #put strains in rows in order of expression
    methyl.pos <- as.numeric(colnames(methyl.mat))
    #plot.methyl.mat(methyl.mat)
        
    if(plot.results){
        if(save.images){
            jpeg(file.path(path, paste0(plot.label, ".Methylation.", gene.name, ".jpg")),
            width = 10, height = 7, units = "in", res = 300)
        }else{
            quartz(width = 10, height = 7)
        }
        
        layout(matrix(c(5,4,1:3), ncol = 1), heights = c(0.5, 0.2, 1,0.2, 1))
        #layout.show(5)
        par(mar = c(0,4,0,4))
        plot.methyl.mat(methyl.mat)

        #add the gene
        par(mar = c(0,4,0,4))
        plot.new()
        plot.window(xlim = c(min(methyl.pos), max(methyl.pos)), ylim = c(0,1))
        if(gene.strand == 1){
            arrows(x0 = gene.start, x1 = gene.end, y0 = 0.5, lwd = 3)
        }else{
            arrows(x0 = gene.end, x1 = gene.start, y0 = 0.5, lwd = 3)
        }
        axis(1)

        #add exons
        par(mar = c(0,4,4,4))
        plot.new()
        plot.window(xlim = c(min(methyl.pos), max(methyl.pos)), 
        ylim = c(0, length(exon.start)))
        for(i in 1:length(exon.start)){
            draw.rectangle(exon.start[i], exon.end[i], i-0.2, i+0.2, fill = "gray")
        }
        mtext("Exons", side = 2)
        mtext(gene.name, side = 3, outer = TRUE, line = -2.5)
    }
    
    island.methyl <- quantify.islands.pos(methyl.mat, plot.bins = FALSE, 
    gene.frac.gap = gene.frac.gap)

    island.bins <- island.methyl[[1]]
    island.values <- island.methyl[[2]]
    
    island.var <- sapply(island.values, var)
    
    if(plot.results){
        par(mar = c(0,4,2,4))
        num.islands <- max(island.bins)
        plot.new()
        plot.window(xlim = c(min(methyl.pos), max(methyl.pos)), ylim = c(0,1))
        for(i in 1:num.islands){
            island.locale <- which(island.bins == i)
            island.pos <- methyl.pos[island.locale]
            min.pos <- min(island.pos)
            max.pos <- max(island.pos)
            draw.rectangle(min.pos, max.pos, 0.25, 0.75, fill = "gray", border.col = "black")
        }
        par(xpd = FALSE)
        mtext("Bins", side = 2)
        }

    get.methyl.cor <- function(expr.mat, methyl.vals, min.rep){
        num.rep <- length(which(!is.na(methyl.vals[1,])))
        if(num.rep < min.rep){
            return(NA)
        }else{
            cor.val <- coef(lm(as.vector(expr.mat)~as.vector(methyl.vals)))[2]
            return(cor.val)
        }
    }

    #methyl.cor <- suppressWaron Expressionnings(sapply(island.values, function(x) cor(x, mean.expr[expr.order], use = "complete")))
    dummy.methyl <- lapply(island.values, 
    function(x) matrix(rep(x, nrow(gene.expr)), nrow = nrow(gene.expr), byrow = TRUE))
    
    methyl.cor <- sapply(dummy.methyl, 
    function(x) get.methyl.cor(gene.expr[,expr.order], x, min.representation))
    names(methyl.cor) <- NULL
    
    if(plot.results){
        plot.new()
        plot.window(xlim = c(min(methyl.pos), max(methyl.pos)), ylim = c(-1, 1))
        for(i in 1:length(methyl.cor)){
            island.locale <- which(island.bins == i)
            island.pos <- methyl.pos[island.locale]
            min.pos <- min(island.pos)
            max.pos <- max(island.pos)
            x <- mean(c(min.pos, max.pos))
            y <- methyl.cor[i]
            segments(x0 = x, y0 = y, y1 = 0, lwd = 3)
        }
        axis(2)
        abline(h = 0)
        mtext("Effect (beta)", side = 2, line = 2.5) 
        if(save.images){dev.off()}
    }

    if(!sliding.window){strains <- strains[1]}

    if(length(strains) > 1){
        if(plot.results){
            if(save.images){
                jpeg(file.path(path, paste0(plot.label, ".Windows.", gene.name, ".jpg")),
                width = 10, height = 5, units = "in", res = 300)
            }else{
                quartz(width = 10, height = 5)
            }
        }
        window.cor <- plot.best.methyl.cor(gene.name, methyl.mat, 
        ordered.expr = gene.expr[,expr.order], gene.start = gene.start, 
        gene.end = gene.end, strand = gene.strand, 
        min.window = min.window, max.window = max.window, plot.results = plot.results)
        
        if(plot.results && save.images){dev.off()}
    }else{
        window.cor <- NA
    }
    
    final.results <- list("gene.info" = gene.info.table, "expr" = gene.expr, 
    "methylation.matrix" = methyl.mat, "island.bins" = island.bins, 
    "island.avg.methylation" = island.values, "island.expr.correlation" = methyl.cor, 
    "window.cor" = window.cor)
    invisible(final.results)

}