#This function plots the correlation between methylation and 
#inbred expression overlayed on the correlation between chromatin
#state and inbred expression

plot.chrom.and.methyl.cor <- function(chrom.cor, methyl.cor, state.colors, plot.label){
	
	
	methyl.start <- unlist(lapply(strsplit(rownames(methyl.cor), ":"), function(x) x[1]))
	methyl.stop <- unlist(lapply(strsplit(rownames(methyl.cor), ":"), function(x) x[2]))	
	methyl.med <- round(rowMeans(cbind(as.numeric(methyl.start), as.numeric(methyl.stop))))
	
	chrom.start <- unlist(lapply(strsplit(rownames(chrom.cor), ":"), function(x) x[1]))
	chrom.stop <- unlist(lapply(strsplit(rownames(chrom.cor), ":"), function(x) x[2]))	
	chrom.med <- round(rowMeans(cbind(as.numeric(chrom.start), as.numeric(chrom.stop))))
	
	xmin <- min(c(methyl.med, chrom.med))
	xmax <- max(c(methyl.med, chrom.med))	
	
	
	layout(matrix(c(1,2), ncol = 1))
	par(mar = c(0,4,4,4))
	plot.new()
	plot.window(ylim = c(-1, 1), xlim = c(xmin, xmax))
	points(x = methyl.med, y = methyl.cor[,1], type = "l", lwd = 2)
	axis(2);mtext("Methylation Correlation", side = 2, line = 2.5)
	abline(h = 0)
	par(mar = c(4,4,0,4))
	plot.new()
	plot.window(ylim = c(-1, 1), xlim = c(xmin, xmax))
	for(i in 1:ncol(chrom.cor)){
		points(chrom.med, chrom.cor[,i], col = state.colors[i], lwd = 2, type = "l")
		}
	axis(2);mtext("Histone State Correlation", side = 2, line = 2.5)
	abline(h = 0)	
	mtext(plot.label, side = 3, outer = TRUE, line = -2)
	
	all.seq <- xmin:xmax
	full.mat <- matrix(NA, nrow = length(all.seq), ncol = (ncol(chrom.cor)+1))
	rownames(full.mat) <- all.seq
	full.mat[as.character(methyl.med),1] <- methyl.cor[,1]
	full.mat[as.character(chrom.med),2:ncol(full.mat)] <- chrom.cor
	all.na <- apply(full.mat, 1, function(x) all(is.na(x)))
	full.mat <- full.mat[which(!all.na),]
	interp <- apply(full.mat, 2, function(x) approx(as.numeric(rownames(full.mat)), x, xout = as.numeric(rownames(full.mat))))
	new.mat <- Reduce(cbind, lapply(interp, function(x) x$y))
	
	# quartz()
	# layout(matrix(c(1,2), ncol = 1))
	# par(mar = c(0,4,4,4))
	# plot.new()
	# plot.window(ylim = c(-1, 1), xlim = c(xmin, xmax))
	# points(x = as.numeric(rownames(full.mat)), y = new.mat[,1], type = "l", lwd = 2)
	# axis(2);mtext("Methylation Correlation", side = 2, line = 2.5)
	# abline(h = 0)
	# par(mar = c(4,4,0,4))
	# plot.new()
	# plot.window(ylim = c(-1, 1), xlim = c(xmin, xmax))
	# for(i in 2:ncol(new.mat)){
		# points(as.numeric(rownames(full.mat)), new.mat[,i], col = state.colors[(i-1)], lwd = 2, type = "l")
		# }
	# axis(2);mtext("Histone State Correlation", side = 2, line = 2.5)
	# abline(h = 0)	
	# mtext(plot.label, side = 3, outer = TRUE, line = -2)

	# all.cor <- apply(new.mat[,2:ncol(new.mat)], 2, function(x) cor(new.mat[,1], x, use = "complete"))
	# plot(new.mat[,1], new.mat[,7])
}