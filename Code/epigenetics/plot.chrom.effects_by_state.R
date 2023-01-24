plot.chrom.effects_by_state <- function(gene.name, expr, covar, chrom.mats, transcript.info, 
transcript.haplotypes, chrom.states, strain.key){
    
    #use the chromatin state object to determine the number of states
    num.states = max(sapply(chrom.mats, 
    function(x) if(length(x) > 1){max(x, na.rm = TRUE)}else{NA}), na.rm = TRUE)
    data(CCcolors)
    chrom.colors <- colors.from.values(1:num.states, use.pheatmap.colors = TRUE)

    gene.id <- transcript.info[which(transcript.info[,"external_gene_name"] == gene.name),"ensembl_gene_id"][1]
    expr.locale <- which(colnames(expr) == gene.id)

    if(length(expr.locale) == 0){
        plot.text("Not in this data set.")
        return(NULL)
    }

    chrom.locale <- which(names(chrom.mats) == gene.id)

    #chromatin scan
    chrom.geno <- get_one_geno(gene.name, transcript.info, transcript.haplotypes, 
    chrom.states, geno_type = "chromatin", strain.key = strain.key)
    
    if(is.na(chrom.geno)){
        plot.text("Not in this data set.")
        return(NULL)
    }

    #check for variation
    test.var <- apply(chrom.geno[[1]], 3, function(x) apply(x, 2, function(y) var(y)))
    no.var <- all(apply(test.var, 1, function(x) length(unique(x))) == 1)

    if(no.var){
        plot.text("No variation in chromatin for this gene.")
        return(NULL)
    }

    #test each chromatin state and each haplotype separately so we 
    #can compare the
    adj.expr <- adjust(expr[,expr.locale,drop=FALSE], covar)
    
    all.chrom.models <- apply(chrom.geno[[1]], 3, function(x) apply(x, 2, function(y) lm(adj.expr~y)))
    all.chrom.r2 <- sapply(all.chrom.models, function(x) sapply(x, function(y) summary(y)$adj.r.squared))
    #chrom.r2.max <- apply(all.r2, 2, max) #maximum R2 by position
    #chrom.r2.max <- apply(all.r2, 1, max) #maximum R2 by state
    all.chrom.coef <- sapply(all.chrom.models, function(x) sapply(x, function(y) coef(y)[2]))

    #get haplotype R2
    geno.geno <- get_one_geno(gene.name, transcript.info, transcript.haplotypes, 
        chrom.states, geno_type = "genotype", strain.key = strain.key)
    all.geno.models <- apply(geno.geno[[1]], 3, function(x) apply(x, 2, function(y) lm(adj.expr~y)))
    all.geno.r2 <- sapply(all.geno.models, function(x) sapply(x, function(y) summary(y)$adj.r.squared))
    all.geno.coef <- sapply(all.geno.models, function(x) sapply(x, function(y) coef(y)[2]))

    #1: chromatin state legend
    #2: chromatin state coefficients
    #3: haplotype coefficient bar plot
    #4: gene coordinates
    #5: chromatin R2 scores
    #6: chromatin state matrix
    layout.mat <- matrix(c(2,1,3,1,5,7,6,4), ncol = 2, byrow = TRUE)
    layout(layout.mat, widths = c(1,0.4), heights = c(1, 0.3, 1, 1))
    #layout.show(7)
    
    #plot state legend, plot 1
    par(mar = c(0,0,4,2))
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,1))
    legend(0,1, fill = chrom.colors, legend = paste("State", 1:num.states))
    
    #plot state coefficients, plots 2 and 3
    par(mar = c(2,4,2,4))
    plot.coef_by_state(all.chrom.coef, gene.id, transcript.info, "chromatin")

    coef.order <- order(all.geno.coef, decreasing = FALSE)
    ordered.coef <- all.geno.coef[coef.order,1]
    ordered.coef.names <- strain.key[match.order(gsub(".y", "", names(ordered.coef)), 
    strain.key[,1], strain.key),1]
    
    #plot DO haplotype coefficients, plot 4
    par(mar = c(1,4,1,4))
    barplot(ordered.coef, col = CCcolors[coef.order], 
    xlab = "", names = ordered.coef.names, las = 2,
    main = "", horiz = TRUE, axes = FALSE)
    axis(3)
    par(xpd = NA)
    mtext("Haplotype Coefficients", side = 3, line = 2.5)
    par(xpd = TRUE)
    
	
    #plot chromatin R2 scores, plot 5
    par(mar = c(0,4,2,4))
    ymax <- max(c(all.chrom.r2, all.geno.r2), na.rm = TRUE)
    ymin <- min(c(all.chrom.r2, all.geno.r2), na.rm = TRUE)
    coord <- as.numeric(colnames(all.chrom.r2))

    state.cols <- colors.from.values(1:num.states, use.pheatmap.colors = TRUE)
    plot.new()
    plot.window(xlim = c(min(coord), max(coord)), ylim = c(ymin, ymax))
    for(s in 1:nrow(all.chrom.r2)){
        points(coord, all.chrom.r2[s,], type = "l", col = state.cols[s])
    }
    axis(2)
    mtext("Variance Explained", side = 2, line = 2.5)
    mtext("Variance Explained", side = 3)

    segments(x0 = min(coord), x1 = max(coord), y0 = all.geno.r2, col = CCcolors)
	par(xpd = NA)
	max.x <- max(coord)
    text(x = max.x, y = all.geno.r2, names(CCcolors), adj = 0)
	par(xpd = TRUE)
    
    
    #plot chromatin state matrix, plot 6
    par(mar = c(0,4,0,4))
    state.order <- match.order(rev(gsub(".y", "", names(ordered.coef))), 
        colnames(chrom.mats[[chrom.locale]]), strain.key)
    #names(ordered.coef); colnames(chrom.mats[[chrom.locale]])[state.order]
    plot.chrom.mat(t(chrom.mats[[chrom.locale]][,state.order]), num.states = num.states)
    
    #plot legend for genetic eQTL LOD score in LOD score plot, plot 7
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,1))
    par(mar = c(2,0,2,0))
    #segments(x0 = 0, x1 = 0.2, y0 = 0.5, lty = 2, 
    #col = "gray", lwd = 3)
    #text(x = 0.25, y = 0.5, "eQTL LOD Score", adj = 0)

}