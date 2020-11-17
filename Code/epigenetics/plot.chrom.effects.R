plot.chrom.effects <- function(gene.name, expr, covar, chrom.mats, transcript.info, 
transcript.haplotypes, chrom.states, strain.key){
    
    num.states = 8
    data(CCcolors)
    chrom.colors <- colors.from.values(1:8, use.pheatmap.colors = TRUE)

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

    chrom.lod <- scan1(chrom.geno, expr[,expr.locale], addcovar = covar)  
    chrom.coef <- scan1coef(genoprobs = chrom.geno, pheno = expr[,expr.locale], 
    addcovar = covar)
    chrom.models <- apply(chrom.geno[[1]], 3, function(x) lm(expr[,expr.locale]~x+covar))
    chrom.r2 <- sapply(chrom.models, function(x) summary(x)$adj.r.squared)

    #scan eQTL
    geno.geno <- get_one_geno(gene.name, transcript.info, transcript.haplotypes, 
    chrom.states, geno_type = "genotype", strain.key = strain.key)
    geno.lod.results <- scan1(geno.geno, expr[,expr.locale], addcovar = covar)
    geno.r2 <- summary(lm(expr[,expr.locale]~geno.geno[[1]][,,1]+covar))$adj.r.squared
    geno.lod <- geno.lod.results[1,1]
    geno.coef <- scan1coef(genoprobs = geno.geno, pheno = expr[,expr.locale], 
    addcovar = covar)[,LETTERS[1:8]]

    #1: chromatin state legend
    #2: chromatin state coefficients
    #3: haplotype coefficient bar plot
    #4: gene coordinates
    #5: chromatin LOD scores
    #6: chromatin state matrix
    layout.mat <- matrix(c(1,0,2,0,3,0,5,7,6,4), ncol = 2, byrow = TRUE)
    layout(layout.mat, widths = c(1,0.4), heights = c(0.3, 1, 0.3, 1, 1))
    #layout.show(7)
    
    #plot state legend, plot 1
    par(mar = c(0,4,0,4))
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,1))
    legend(0,0.5, fill = chrom.colors, legend = paste("State", 1:num.states), 
    horiz = TRUE)
    
    #plot state coefficients, plots 2 and 3
    par(mar = c(2,4,2,4))
    plot.coef(chrom.coef, gene.id, transcript.info, "chromatin")

    ordered.coef <- sort(geno.coef, decreasing = FALSE)
    ordered.coef.names <- strain.key[match.order(names(ordered.coef), 
    strain.key[,1], strain.key),1]
    
    #plot DO haplotype coefficients, plot 4
    par(mar = c(1,4,1,4))
    barplot(ordered.coef, col = CCcolors[order(geno.coef, decreasing = FALSE)], 
    xlab = "", names = ordered.coef.names, las = 2,
    main = "", horiz = TRUE, axes = FALSE)
    axis(3)
    par(xpd = NA)
    mtext("Haplotype Coefficients", side = 3, line = 2.5)
    par(xpd = TRUE)
    
	
    #plot chromatin LOD scores, plot 5
    par(mar = c(0,4,2,4))
    ymax <- max(c(chrom.lod[,1], geno.lod), na.rm = TRUE)
    plot(chrom.lod[,1], type = "l", ylab = "LOD score", axes = FALSE,
    xlab = "", ylim = c(0, ymax))
    abline(h = geno.lod, col = "gray", lwd = 3, lty = 2)
    axis(2)

	par(xpd = NA)
	max.x <- nrow(chrom.lod)*1.15
    text(x = max.x, y = geno.lod, "eQTL LOD Score", adj = 0)
	par(xpd = TRUE)
    
    
    #plot chromatin state matrix, plot 6
    par(mar = c(0,4,0,4))
    state.order <- match.order(rev(names(ordered.coef)), 
    colnames(chrom.mats[[chrom.locale]]), strain.key)
    #names(ordered.coef); colnames(chrom.mats[[chrom.locale]])[state.order]
    plot.chrom.mat(t(chrom.mats[[chrom.locale]][,state.order]))
    
    #plot legend for genetic eQTL LOD score in LOD score plot, plot 7
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,1))
    par(mar = c(2,0,2,0))
    #segments(x0 = 0, x1 = 0.2, y0 = 0.5, lty = 2, 
    #col = "gray", lwd = 3)
    #text(x = 0.25, y = 0.5, "eQTL LOD Score", adj = 0)

}