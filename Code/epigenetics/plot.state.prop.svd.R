plot.state.prop.svd <- function(state.prop.mat, gene.name, gene.info.table, group.mean.expr, state.key){

	state.colors <- c("blue", "green", "gray", "yellow", "red", "purple")
	state.prop.mat <- apply(state.prop.mat, 2, function(x) x - mean(x))
	state.prop.eig <- svd(state.prop.mat)

	layout(matrix(c(1,2), ncol = 1), heights = c(1,2))
	par(mar = c(1,4,2,2))
	singular.vals <- state.prop.eig$d
	eigen.vals <- singular.vals^2
	var.accounted <- eigen.vals/sum(eigen.vals)
	barplot(var.accounted, main = "Variance Accounted")
	par(mar = c(4,4,0,2))
	plot(state.prop.eig$v[,1], state.prop.eig$v[,2], col = state.colors, pch = 16, xlab = "PC1", ylab = "PC2")
	
	gene.locale <- which(gene.info.table[,2] == gene.name)
	gene.id <- gene.info.table[gene.locale,1]
	id.locale <- which(names(group.mean.expr) == gene.id)
	gene.expr <- group.mean.expr[[id.locale]]
	
	plot(state.prop.eig$v[,1], gene.expr)
	
	
	invisible(state.prop.eig)
	
	}