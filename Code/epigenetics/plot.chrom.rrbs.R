#This function plots a chromatin matrix and a methylation matrix
#on the same coordinates
#I have added an optional argument for sequence information as well
#if gene.vars is supplied, and a strain order other than the default
#is provided, col.table must also be supplied.

plot.chrom.rrbs <- function(methyl.mat, chrom.mat, gene.name, gene.info, 
gene.frac.gap = 0.01, gene.vars = NULL, col.table = NULL){

    gene.locale <- which(gene.info[,"external_gene_name"] == gene.name)
    if(length(gene.locale) == 0){stop("can't find ", gene.name)}

    methyl.pos <- as.numeric(rownames(methyl.mat))
    chrom.pos <- as.numeric(rownames(chrom.mat))
    all.pos <- sort(unique(c(methyl.pos, chrom.pos)))

    gene.chr <- unique(gene.info[gene.locale,"chromosome_name"])
    gene.start <- unique(gene.info[gene.locale,"start_position"])
    gene.end <- unique(gene.info[gene.locale,"end_position"])
    gene.strand <- unique(gene.info[gene.locale,"strand"])
    exon.start <- unique(gene.info[gene.locale,"exon_chrom_start"])
    exon.end <- unique(gene.info[gene.locale,"exon_chrom_end"])

    gene.info.table <- list("gene.name" = gene.name, "Chr" = gene.chr, 
    "gene_start" = gene.start, "gene_end" = gene.end, "gene_strand" = gene.strand, 
    "exon_start" = exon.start, "exon_end" = exon.end)

    if(is.null(gene.vars)){
        layout(matrix(c(1:3), ncol = 1), heights = c(1, 0.2, 1))
    }else{
        layout(matrix(c(1:4), ncol = 1), heights = c(1, 0.2, 1, 1))
    }

    par(mar = c(0,4,4,4))
    plot.methyl.mat(t(methyl.mat), xlim = c(min(all.pos), max(all.pos)))

     par(mar = c(2,4,0,4))
        plot.new()
        plot.window(xlim = c(min(methyl.pos), max(methyl.pos)), ylim = c(0,1))
        if(gene.strand == 1){
            arrows(x0 = gene.start, x1 = gene.end, y0 = 0.5, lwd = 3)
        }else{
            arrows(x0 = gene.end, x1 = gene.start, y0 = 0.5, lwd = 3)
        }
        axis(1)

        #add exons
       # par(mar = c(0,4,4,4))
       # plot.new()
       # plot.window(xlim = c(min(methyl.pos), max(methyl.pos)), 
       # ylim = c(0, length(exon.start)))
       # for(i in 1:length(exon.start)){
       #     draw.rectangle(exon.start[i], exon.end[i], i-0.2, i+0.2, fill = "gray")
       # }
       # mtext("Exons", side = 2)
        mtext(gene.name, side = 3, outer = TRUE, line = -2.5)
    
        islands <- quantify.islands.pos(t(methyl.mat), gene.frac.gap = gene.frac.gap)
        island.pos <- cbind(islands$island.position, as.numeric(rownames(methyl.mat)))

        par(mar = c(0,4,0,4))
        rev.mat <- chrom.mat[,ncol(chrom.mat):1]
        plot.chrom.mat(state.mat = t(rev.mat), num.states = 8, 
        island.bins = island.pos, xlim = c(min(all.pos), max(all.pos)))

        if(!is.null(gene.var)){
            par(mar = c(4,4,0,4))
            plot.variants(gene.vars, xlim = c(min(all.pos), max(all.pos)),
            strain.order = colnames(rev.mat), col.table = col.table)        
            plot.island.obj(islands, add = TRUE)
        }

    #

}