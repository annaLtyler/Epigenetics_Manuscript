# chr = 6; start.coord = 121650075;stop.coord = 121654314

get.nearest.gene <- function(gene.table, chr, start.coord, stop.coord, upstream.buffer = 5000, downstream.buffer = 5000){
	
	#make sure the table only includes protein coding genes
	all.gene.locale <- which(gene.table[,2] == "protein coding gene")
	subtable <- gene.table[all.gene.locale,]
	
	chr.locale <- which(subtable[,5] == chr)
	chr.subtable <- subtable[chr.locale,]
	gene.starts <- chr.subtable[,6] - upstream.buffer
	gene.ends <- chr.subtable[,7] + downstream.buffer
		
	#find everywhere the mark start position is within the gene boundary
	mark.start.inc <- intersect(which(gene.starts <= stop.coord), which(gene.ends >= start.coord))
	#find everywhere the mark end position is within the gene boundary
	mark.stop.inc <- intersect(which(gene.starts <= start.coord), which(gene.ends >= start.coord))
	
	#overlaps are where either the mark start or end are within the gene boundary
	overlaps <- unique(c(mark.start.inc, mark.stop.inc))

	overlap.genes <- unique(chr.subtable[overlaps,3])

	return(overlap.genes)

	# overlap.starts <- subtable[overlaps,6]
	# overlap.ends <- subtable[overlaps,7]

	# plot.min <- min(c(overlap.starts, overlap.ends))
	# plot.max <- max(c(overlap.starts, overlap.ends))
	# plot.width = plot.max - plot.min

	# plot.new()
	# plot.window(xlim = c(plot.min, ((plot.max + (plot.width*0.06)))), ylim = c(0,length(overlap.genes)))
	# par(xpd = TRUE)
	# segments(start.coord, 0, end.coord, 0, col = "black", lwd = 3)
	# text(x = (plot.max + (plot.width*0.03)), y = 0, "mark", adj = 0)
	# for(i in 1:length(overlap.genes)){
		# segments(overlap.starts[i], i, overlap.ends[i], i, col = "blue", lwd = 3)
		# text(x = (plot.max + (plot.width*0.03)), y = i, overlap.genes[i], adj = 0)
		# }
	# par(xpd = FALSE)
}