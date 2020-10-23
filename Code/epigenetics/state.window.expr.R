#This function correlates every point along a gene's state matrix with expression
#gene.exon.table is required for biologically informed binning of genes

state.window.expr <- function(gene.name, gene.info.table, gene.exon.table = NULL, 
chrom.mats, group.mean.expr, state.key, window.gene.prop = 0.1, gap.prop = 0.005, 
bin.by = c("window", "bio"), collapse.exons = TRUE, plot.results = TRUE){

	state.col <- colors.from.values(1:num.states, use.pheatmap.colors = TRUE, 
	global.color.scale = TRUE, global.min = 1, global.max = num.states)

	#====================================================
	# internal functions
	#====================================================
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

	bin.by <- bin.by[1]
	
	gene.locale <- which(gene.info.table[,"external_gene_name"] == gene.name)[1]
	if(length(gene.locale) == 0){
		# cat("I couldn't find ", gene.name, "\n")
		return(matrix(NA, ncol = 6, nrow = 1))
		}
		
	# if(length(gene.locale) > 1){stop("I found more than one gene with this ID")}

	for(g in 1:length(gene.locale)){

	tss <- gene.info.table[gene.locale[g],"start_position"]	
	tes <- gene.info.table[gene.locale[g],"end_position"]	
	strand <- gene.info.table[gene.locale[g],6]
	
	gene.id <- gene.info.table[gene.locale[g],1]
	
	state.locale <- which(names(chrom.mats) == gene.id)
	expr.locale <- which(names(group.mean.expr) == gene.id)
	
	if(length(expr.locale) == 0){
		# cat("I couldn't find expression values for ", gene.name, "\n")
		return(paste("no expression", gene.name))
		}

	if(length(state.locale) == 0){
		# cat("I couldn't find chromatin states for ", gene.name, "\n")
		return(paste("no chromatin", gene.name))
		}
	
	expr <- group.mean.expr[[expr.locale]]
	states <- chrom.mats[[state.locale]]
	
	if(all(is.na(expr)) || all(is.na(states))){
		# cat("No expression data for", gene.name, "\n")
		return(paste("no expression", gene.name))
		}
	
	max.states <- nrow(state.key)
	
	#use a sliding window across the states and calculate
	#the correlation between the proportion of each state
	#in each window and the mean inbred expression
	bin.cor <- function(state.mat, bins.which, expr.vals){
		
		if(is.null(bins.which)){
			return(rep(NA, max.states))
			}
		sub.mat <- state.mat[,bins.which,drop=FALSE]

		#calculate proportions of each state in the bin
		state.props <- lapply(1:max.states, function(st) apply(sub.mat, 1, function(x) length(which(x == st))/length(x)))
		state.prop.mat <- Reduce("rbind", state.props)
		
		state.cor <- apply(state.prop.mat, 1, function(x) cor(x, expr.vals))
		names(state.cor) <- 1:max.states
		return(state.cor)
		}
	

	# pdf("Diff.windows.pdf", width = 10, height = 5)
	# for(window.size in (nrow(states)/4):(nrow(states)-1)){	
	# report.progress(window.size, (nrow(states)-1))

	if(bin.by == "window"){
		window.size <- round(nrow(states)*window.gene.prop)
		gap.size <- round(nrow(states)*gap.prop)
		if(gap.size < 1){gap.size = 1}
		bins <- sliding.window.el(1:nrow(states), window.size, gap.size)
		bin.names <- unlist(lapply(bins, function(x) paste0(rownames(states)[x[1]], ":", rownames(states)[tail(x, 1)])))
		}else{
		gene.exon.locale <- which(gene.exon.table[,2] == gene.name)
		exon.table <- gene.exon.table[gene.exon.locale,,drop=FALSE]
		tss <- gene.info.table[gene.locale[g],"start_position"]
		gene.end <- gene.info.table[gene.locale[g],"end_position"]		
		bins <- bin.gene(states, strand, tss, gene.end, exon.table, collapse.exons)	
		bin.pos <- lapply(bins, function(x) paste0(rownames(states)[x[1]], ":", rownames(states)[tail(x, 1)]))
		bin.names <- names(bins)
		}
	
	all.state.cor <- lapply(bins, function(x) bin.cor(t(states), x, expr))
	all.state.cor.mat <- Reduce("rbind", all.state.cor)
	rownames(all.state.cor.mat) <- bin.names
	
	
	if(bin.by == "window"){

	if(strand == -1){all.state.cor.mat <- all.state.cor.mat[nrow(all.state.cor.mat):1,]}

	
	
	if(plot.results){
		layout(matrix(c(1:4), ncol = 2, byrow = TRUE), widths = c(1,0.3))
		par(mar = c(3,4,4,0))
		plot.new()
		plot.window(ylim = c(-1, 1), xlim = c(0, nrow(all.state.cor.mat)))
		for(i in 1:ncol(all.state.cor.mat)){
			points(1:nrow(all.state.cor.mat), all.state.cor.mat[,i], type = "l", col = state.colors[i], lwd = 3)
			}
		axis(2)
		mtext("Correlation", side = 2, line = 2.5)
		
		tss.bins <- get.in.bin(rownames(all.state.cor.mat), tss, strand)
		tes.bins <- get.in.bin(rownames(all.state.cor.mat), tes, strand)
		par(xpd = TRUE)
		arrows(min(tss.bins), -1.1, max(tes.bins), -1.1, length = 0.1, lwd = 3)
		
		plot.new()
		plot.window(xlim = c(0,1), ylim = c(0,1))
		par(mar = c(0,0,0,0))
		legend(-0.6, 0.5, col = state.colors, lty = 1, legend = state.key[,2])	
		# }
		# dev.off()
		}
	
	best.of.each.state <- apply(all.state.cor.mat, 2, function(x) x[which(abs(x) == max(abs(x), na.rm = TRUE))[1]])
	
	best.position <- apply(all.state.cor.mat, 2, function(x) which(abs(x) == max(abs(x), na.rm = TRUE))[1])
	
	max.window <- lapply(best.position, function(x) bins[[x]]/nrow(states))

	# par(mar = c(4,4,2,2))
	# expr.order <- match(names(expr), colnames(states))
	# col.names <- unlist(lapply(strsplit(colnames(states), "/"), function(x) x[1]))
	# mat.to.plot <- states[nrow(states):1,expr.order]
	# min.split <- 0.5
	# max.split <- 6.5
	# imageWithText(mat = t(states), show.text = FALSE, row.names = col.names[expr.order], col.text.rotation = 0, main = "Methylation States in Controls", split.at.vals = TRUE, split.points = c(min.split:max.split), col.scale = c("blue", "green", "gray", "yellow", "red", "purple"), light.dark = "d")

	if(plot.results){
		par(mar = c(3,4,4,0))
		plot.new()
		plot.window(xlim = c(0,1), ylim = c(-1,1))
		points(x = c(0, 1), y = c(0,0), type = "l", col = "darkgray", lwd = 3, lty = 2)
		for(i in 1:length(best.position)){
			points(x = max.window[[i]], y = rep(best.of.each.state[[i]], length(max.window[[i]])), col = state.colors[i], type = "l", lwd = 3)
			}
			axis(1);axis(2)
	
		plot.new()
		plot.window(xlim = c(0,1), ylim = c(0,1))
		par(mar = c(0,0,0,0))
		legend(-0.6, 0.5, col = state.colors, lty = 1, legend = state.key[,2])	
	
		# dev.off()
		mtext(gene.name, side = 3, outer = TRUE, line = -2.5)
		}
		
	# return(rbind(best.of.each.state, best.position/nrow(all.state.cor.mat)))
	#end code for binning by window
	}else{
		#begin code for biological binning
		upstream.locale <- which(rownames(all.state.cor.mat) == 'upstream')
		downstream.locale <- which(rownames(all.state.cor.mat) == 'downstream')
		promoter.locale <- which(rownames(all.state.cor.mat) == 'promoter')
		intron.locale <- grep("intron", rownames(all.state.cor.mat))
		exon.locale <- grep("exon", rownames(all.state.cor.mat))
		
		
		if(plot.results){
		layout(matrix(c(1:6), ncol = 2, byrow = TRUE))
		par(mar = c(2,2,2,2))
		barplot(all.state.cor.mat[upstream.locale,], col = state.colors, ylim = c(-1,1), names = NA, main = "Upstream Region"); abline(h = 0)
		

		barplot(all.state.cor.mat[promoter.locale,], col = state.colors, ylim = c(-1,1), names = NA, main = "Promoter Region"); abline(h = 0)

		
		if(collapse.exons){
			barplot(all.state.cor.mat[intron.locale,], col = state.colors, ylim = c(-1,1), names = NA, main = "Introns"); abline(h = 0)
				}else{
			plot.new()
			plot.window(xlim = c(1,max.states), ylim = c(-1, 1))	
			intron.mat <- all.state.cor.mat[intron.locale,]
			for(i in 1:ncol(intron.mat)){
				points(x = rep(i, nrow(intron.mat)), intron.mat[,i], col = state.colors[i], pch = 16)
				}
			axis(2)
			abline(h = 0)
			# par(xpd = TRUE)
			# text(y = rep(-1.2, max.states), x = 1:max.states, labels = state.key[,2])
			mtext("Introns", side = 3, outer = FALSE)
			# par(xpd = FALSE)
			}
		}
		
		if(plot.results){
			if(collapse.exons){
			
			barplot(all.state.cor.mat[exon.locale,], col = state.colors, ylim = c(-1,1), names = NA, main = "Exons"); abline(h = 0)
	
			}else{
			plot.new()
			plot.window(xlim = c(1,max.states), ylim = c(-1, 1))	
			exon.mat <- all.state.cor.mat[exon.locale,]		
			for(i in 1:ncol(exon.mat)){
				points(x = rep(i, nrow(exon.mat)), exon.mat[,i], col = state.colors[i], pch = 16)
				}
			axis(2)
			abline(h = 0)
			# par(xpd = TRUE)
			# text(y = rep(-1.2, max.states), x = 1:max.states, labels = state.key[,2])
			mtext("Exons", side = 3, outer = FALSE)
			par(xpd = FALSE)
			}
			
			
			barplot(all.state.cor.mat[downstream.locale,], col = state.colors, ylim = c(-1,1), main = "Downstream Region", names = NA); abline(h = 0)
	
			plot.new()
			plot.window(xlim = c(0, 1), ylim = c(0,1))
			legend(0.25, 0.8, legend = state.key[,2], fill = state.colors, cex = 1.5)
	
			mtext(gene.name, outer = TRUE, line = -1.5)
			}
		}
	} #end looping through multiple gene ids for one gene name
	return(all.state.cor.mat)
	
}