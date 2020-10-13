dist.to.nearest.TSS <- function(gene.info.table, start.pos, end.pos, chr){
	
	chr.locale <- which(gene.info.table[,3] == chr)
	chr.subtable <- gene.info.table[chr.locale,]
	start.to.start = get.nearest.pt(chr.subtable[,4], start.pos)
	start.to.end = get.nearest.pt(chr.subtable[,4], end.pos)	
	
	nearest.pts <- unique(c(start.to.start, start.to.end))
	
	gene.table <- chr.subtable[nearest.pts,,drop=FALSE]
	
	dist.to.state.start <- unlist(lapply(gene.table[,4], function(x) x - start.pos))
	dist.to.state.end <- unlist(lapply(gene.table[,4], function(x) x - end.pos))
	
	dist.to.state <- apply(cbind(dist.to.state.start, dist.to.state.end), 1, function(x) x[which.min(abs(x))])

	final.table <- cbind(gene.table, dist.to.state)
	return(final.table)
}