expr.prop.state.cor <- function(chrom.state.props, group.mean.expr, name.color.key){
	
	shared.genes <- unique(c(names(chrom.state.props), names(group.mean.expr)))
	num.states <- nrow(chrom.state.props[[1]])
	num.strains <- ncol(chrom.state.props[[1]])
	
	chrom.strains <- name.color.key[match(colnames(chrom.state.props[[1]]), name.color.key[,4]),1]
	expr.strains <- name.color.key[match(names(group.mean.expr[[1]]), name.color.key[,4]),1]
	shared.strains <- intersect(chrom.strains, expr.strains)
	chrom.strain.order <- match(shared.strains, chrom.strains)
	expr.strain.order <- match(shared.strains, expr.strains)
	
	gene.cor.mat <- matrix(NA, nrow = length(shared.genes), ncol = num.states)
	gene.p.mat <- matrix(NA, nrow = length(shared.genes), ncol = num.states)
	rownames(gene.cor.mat) <- rownames(gene.p.mat) <- shared.genes
	colnames(gene.cor.mat) <- colnames(gene.p.mat) <- 1:num.states

	for(i in 1:length(shared.genes)){
		report.progress(i, length(shared.genes))
		gene.states.locale <- which(names(chrom.state.props) == shared.genes[i])
		gene.expr.locale <- which(names(group.mean.expr) == shared.genes[i])
		if(!is.na(sum(chrom.state.props[[gene.states.locale]])) && sum(chrom.state.props[[gene.states.locale]]) > 0){
			result <- apply(chrom.state.props[[gene.states.locale]], 1, function(x) cor.test(x[chrom.strain.order], group.mean.expr[[gene.expr.locale]][expr.strain.order]))
			gene.cor.mat[i,] <- unlist(lapply(result, function(x) x$estimate))
			gene.p.mat[i,] <- unlist(lapply(result, function(x) x$p.value))
			}
		}
		
		final.result <- list(gene.cor.mat, gene.p.mat)
		names(final.result) <- c("Correlations.by.State", "P.values.by.States")
		return(final.result)
	
}