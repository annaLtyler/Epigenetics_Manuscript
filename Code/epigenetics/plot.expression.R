#This function plots inbred gene expression

plot.expression <- function(gene.name, gene.info.table, rna.seq, C.or.T = "C", col.table){
	
	gene.locale <- which(gene.info.table[,2] == gene.name)
	gene.id <- gene.info.table[gene.locale,1]
	id.locale <- which(names(group.mean.expr) == gene.id)
	
	gene.expr <- get.gene.expr2(gene.name, gene.info.table, rna.seq, C.or.T)	
	cols <- col.table[match(names(gene.expr), col.table[,4]),3]
	expr.mean <- unlist(lapply(gene.expr, function(x) mean(x[which(is.finite(x))])))
	expr.order <- order(expr.mean)
	
	ymin <- min(unlist(gene.expr)[which(is.finite(unlist(gene.expr)))]);ymax <- max(unlist(gene.expr)[which(is.finite(unlist(gene.expr)))])

	plot.new()
	plot.window(xlim = c(1,length(expr.mean)), ylim = c(ymin, ymax))
	for(i in 1:length(gene.expr)){
		xval <- rep(i, length(gene.expr[[i]]))
		yval <- gene.expr[[expr.order[i]]]
		points(x = xval, y = yval, col = cols[expr.order[i]], pch = 16)
		}

	segments(x0 = 1:length(expr.mean)-0.2, y0 = expr.mean[expr.order], x1 = 1:length(expr.mean)+0.2, y1 = expr.mean[expr.order], lwd = 3, col = cols[expr.order])

	axis(2)
	plot.height <- ymax - ymin
	par(xpd = TRUE)
	text(x = 1:length(gene.expr), y = rep(ymin-(plot.height*0.1), length(expr.mean)), labels = names(gene.expr)[expr.order], col = cols[expr.order], font = 2)
	par(xpd = FALSE)
	mtext(gene.name)

	return(gene.expr)
	
}