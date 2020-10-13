#This function looks at correlation between DNA methylation in a sliding 
#window and inbred expression
#gene.name = "Corin"; gene.info.table = rnaseq.gene.info; upstream.buffer = 1000; downstream.buffer = 1000;window.gene.prop = 0.1; gap.prop = 0.05


methylation_expression_window <- function(smoothed.methyl, gene.expr, gene.name, col.table){
	
	mean.expr <- unlist(lapply(gene.expr, function(x) mean(x[which(is.finite(x))], na.rm = TRUE)))
	
	mean.order <- order(mean.expr)
	
	methyl.order <- match(col.table[,1], colnames(smoothed.methyl))
	expr.order <- match(col.table[,4], names(mean.expr))
	# cbind(colnames(smoothed.methyl)[methyl.order], names(mean.expr)[expr.order])

	ordered.methyl <- smoothed.methyl[,methyl.order[mean.order]]
	ordered.expr <- mean.expr[expr.order[mean.order]]

	all.state.cor <- apply(ordered.methyl, 1, function(x) cor(x, ordered.expr))
	
	plot(all.state.cor, type = "l", ylim = c(-1, 1), main = gene.name); abline(h = 0)
	mtext("Correlation between smoothed methylation coverage and inbred expression", side = 3)

	result <- matrix(all.state.cor, ncol = 1)
	rownames(result) <- rownames(smoothed.methyl)
	return(result)
}