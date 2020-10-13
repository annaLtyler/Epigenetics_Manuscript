#This function plots a smoothed methylation table

plot.smooth.methyl.table <- function(methyl.table, gene.name, calc.type = c("coverage", "mean"), window.gene.prop = 0.1, gap.prop = 0.005, col.table, gene.info.table){
	
		calc.type <- calc.type[1]
	
		window.size <- round(nrow(methyl.table)*window.gene.prop)
		gap.size <- round(nrow(methyl.table)*gap.prop)
		if(gap.size < 1){gap.size = 1}
		bins <- sliding.window.el(1:nrow(methyl.table), window.size, gap.size)
		bin.names <- unlist(lapply(bins, function(x) paste0(rownames(methyl.table)[x[1]], ":", rownames(methyl.table)[tail(x, 1)])))

		smooth.methyl <- function(methyl.row, bins, calc.coverage = FALSE){
			if(calc.coverage){methyl.row[which(is.na(methyl.row))] <- 0}
			bin.vals <- unlist(lapply(bins, function(x) mean(methyl.row[x], na.rm = TRUE)))
			#plot(bin.vals, type = "l")
			return(bin.vals)
			}
	
		get.in.bin <- function(bin.names, gene.pt, strand){
			split.bins <- strsplit(bin.names, ":")
			bin.start <- unlist(lapply(split.bins, function(x) as.numeric(x[1])))
			bin.end <- unlist(lapply(split.bins, function(x) as.numeric(x[2])))
			if(strand == 1){
				after.start <- which(bin.start <= gene.pt)
				before.end <- which(bin.end >= gene.pt)
				}else{
				after.start <- which(bin.start >= gene.pt)
				before.end <- which(bin.end <= gene.pt)
				}
			bins.containing <- intersect(after.start, before.end)
			return(bins.containing)
			}
	
	
		if(calc.type == "mean"){
		all.smooth <- apply(methyl.table, 2, function(x) smooth.methyl(x, bins, FALSE));ymax = 100;gene.y = -10
		}else{
		all.smooth <- apply(methyl.table, 2, function(x) smooth.methyl(x, bins, TRUE));ymax = max(all.smooth);gene.y <- ymax*-0.05
		}
		
		plot.new()
		plot.window(xlim = c(1,nrow(all.smooth)), ylim = c(0, ymax))	
		for(i in 1:ncol(all.smooth)){
			col <- col.table[match(colnames(methyl.table)[i], col.table[,1]),3]
			points(all.smooth[,i], col = col, lwd = 3, type = "l")
			}
		axis(2)
		
		gene.locale <- which(gene.info.table[,2] == gene.name)
		tss <- gene.info.table[gene.locale,"start_position"]
		tes <- gene.info.table[gene.locale,"end_position"]		
		gene.strand <- gene.info.table[gene.locale,"strand"]		
		bins.with.tss <- get.in.bin(bin.names, tss, gene.strand)
		bins.with.tes <- get.in.bin(bin.names, tes, gene.strand)
		
		par(xpd = TRUE)
		arrows(min(bins.with.tss), gene.y, max(bins.with.tes), gene.y, length = 0.1, lwd = 3)
		par(xpd = FALSE)		
		mtext(gene.name)
		
		rownames(all.smooth) <- bin.names
		invisible(all.smooth)
}