#This function plots stats for genes in a list and 
#puts the results in a pdf with the given label
  	
plot.one.gene <- function(do.eqtl, rnaseq.gene.info, rna.seq, all.bed, col.table, 
gene.name, T.or.C,
start.feature, end.feature, upstream.buffer, downstream.buffer, get.snps = FALSE, 
separate.windows = FALSE, dim.rd = c("mds", "state.prop"), state = NULL, total.states = 8, 
show.state.numbers = FALSE,state.weights = NULL){


	qtl.ids <- names(do.eqtl)
	id <- unique(rnaseq.gene.info[which(rnaseq.gene.info[,"external_gene_name"] == gene.name),1])
	
	if(length(id) == 0){
		plot.new()
		plot.window(xlim = c(0, 1) ,ylim = c(0,1))
		text(0.5, 0.5, paste("There are no eQTL data for", gene.name))	
		return(NULL)
	}

	if(length(id) > 0){
		gene.stats <- get.gene.info(rnaseq.gene.info, rna.seq, all.bed, id, T.or.C,
		start.feature, end.feature, upstream.buffer, downstream.buffer, 
		get.snps = get.snps)
		
		if(is.null(gene.stats)){
			plot.text("No data for gene")
			return(NULL)
		}

		#gene.stats$control.expression <- lapply(gene.stats$control.expression, 
		#function(x) log2(x+1))
		
		plot.expr.and.states(gene.info = gene.stats$gene.info, state.mat = 
		gene.stats$state.matrix, gene.expr = gene.stats$control.expression, 
		snp.tree = gene.stats$snp.tree, chrom.tree = gene.stats$chromatin.state.tree, 
		dim.rd = dim.rd, state = state, num.states = total.states, 
		show.states = show.state.numbers, state.weights = state.weights, 
		col.table = col.table)
		
		gene.start <- gene.stats[[1]][1,4]
		gene.end <- gene.stats[[1]][1,5]
		qtl.locale <- which(qtl.ids == id)	
		if(length(qtl.locale) > 0){
		if(separate.windows){quartz()}
		plot.eqtl(coef.list = do.eqtl[[qtl.locale]][[2]], 
		lod.score = do.eqtl[[qtl.locale]][[3]], label = gene.name, 
		gene.coord = mean(c(gene.start, gene.end)))

		#finally, plot the correlation between the most cis do.eqtl coefficients and the 
		#inbred expression
		inbred.expr <- sapply(gene.stats$control.expression, mean)
		scaled.expr <- scale(inbred.expr)

		chr.coord <- as.numeric(sapply(strsplit(rownames(do.eqtl[[qtl.locale]][[2]]), "_"), 
		function(x) x[2]))
		gene.info = gene.stats$gene.info
		nearest.eqtl <- unique(sapply(gene.info[,"start_position"], function(x) 
		get.nearest.pt(chr.coord, x)))[1]
		cis.coef <- do.eqtl[[qtl.locale]][[2]][nearest.eqtl,]
		coef.order <- order.strains(names(inbred.expr), names(cis.coef), col.table)
		col.order <- order.strains(names(inbred.expr), col.table[,1], col.table)
		
		model <- lm(scaled.expr~cis.coef[coef.order])
		cor.stats <- cor.test(scaled.expr, cis.coef[coef.order], use = "complete")
		r <- signif(cor.stats$estimate, 2)
		p <- signif(cor.stats$p.value, 2)
		
		if(separate.windows){quartz()}
			layout(matrix(c(1, rep(0, 3)), ncol = 2))
			par(mar = c(6,6,4,2))
			plot(cis.coef[coef.order], scaled.expr, pch = 16, col = col.table[col.order,3],
			cex = 2, xlab = "DO cis eQTL Coefficient", ylab = "Inbred Epression", cex.lab = 1.5,
			cex.axis = 1.5, main = paste0("r = ", r, ", p = ", p))
			abline(model)
		}
	}
invisible(gene.stats)
}	
	
