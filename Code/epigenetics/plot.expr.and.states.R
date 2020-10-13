#This function plots multiple high-level objects onto a single
#page, expression and SNP distance trees, expression levels, and
#ChromHMM states
#dim.rd is the dimension reduction method. This can be either 
#multi-dimensional scaling (mds) or using the proportion of a
#state (state.prop). If dim.rd == state.prop, the state must 
#be specified.

plot.expr.and.states <- function(gene.info, state.mat, gene.expr, snp.tree, chrom.tree, 
dim.rd = c("mds", "state.prop"), state = NULL, num.states = 15, show.states = FALSE,
state.weights = NULL, col.table){
		
		dim.rd <- dim.rd[1]
	
		if(dim.rd == "mds"){state <- NULL}

		strains.included <- colnames(state.mat)

		# if(is.null(state)){dim.rd <- "mds"}
		# if(!is.null(state)){dim.rd <- "state.prop"}

		#The first step is to get everything in the same order.
		#otherwise things get very confusing.
		#here we sort expression by mean value. We then order the
		#state matrix and strain colors to match

		
		inc.strain.locale <- order.strains(strains.included, names(gene.expr), col.table)
		gene.expr <- gene.expr[inc.strain.locale]
		expr.order <- order(sapply(gene.expr, mean), decreasing = FALSE)
		ordered.expr <- gene.expr[expr.order]
		#boxplot(ordered.expr)
		
		#reorder the color table to get the colors in order of expression
		strain.order <- order.strains(names(ordered.expr), col.table[,1], col.table)
		ordered.col.table <- col.table[strain.order,] #color.table ordered by expression
		ordered.colors <- ordered.col.table[,3]

		#reorder the states to be in the same order as the expression
		state.order <- order.strains(names(ordered.expr), colnames(state.mat), col.table)
		state.order <- state.order[which(!is.na(state.order))]
		#cbind(colnames(state.mat)[state.order], names(ordered.expr))
		ordered.states <- state.mat[,state.order]		
		
		ordered.names <- ordered.col.table[,2]

		if(is.null(snp.tree) && is.null(chrom.tree)){
			layout.mat <- matrix(c(1,1,2,2,3,4), ncol = 2, byrow = TRUE)
			layout(layout.mat, heights = c(1.5, 1, 1.5))
			# layout.show(4)
		}else{
			layout.mat <- matrix(c(1,2,3,3,4,4), ncol = 2, byrow = TRUE)
			layout(layout.mat, heights = c(1.2, 1, 1))
		}


		if(!is.null(snp.tree) || !is.null(chrom.tree)){
			par(mar = c(1,2,6,0))
			if(!is.null(snp.tree)){
				par(xpd = TRUE)
				plot(snp.tree, main = "SNP comparisons", tip.col = col.table[,3], 
				cex = 1.1, font = 2, type = "unrooted", lab4ut = "axial")
				par(xpd = FALSE)
			}
			if(!is.null(chrom.tree) && length(chrom.tree$tip.label) > 0){
				par(xpd = TRUE)
				plot(chrom.tree, main = "Chromatin comparisons", type = "unrooted", 
				tip.col = col.table[match(chrom.tree$tip.label, col.table[,2]),3], cex = 1.1, 
				font = 2, lab4ut = "axial")
				par(xpd = FALSE)
			}else{
				plot.new()
				plot.window(xlim = c(0,1), ylim = c(0,1))
				text(0.5, 0.5, labels = "No Tree")
			}
		}
		

		par(mar = c(4,4,6,2))
		rownames(ordered.states) <- NULL
		imageWithText(mat = t(ordered.states), show.text = show.states, 
		row.names = colnames(ordered.states), col.text.rotation = 0, 
		main = "Methylation States in Controls", use.pheatmap.colors = TRUE, 
		global.color.scale = TRUE, global.min = 1, global.max = num.states)
		par(xpd = TRUE)

		#draw a legend for the states
		state.col <- colors.from.values(1:num.states, use.pheatmap.colors = TRUE, 
		global.color.scale = TRUE, global.min = 1, global.max = num.states)

		plot.dim <- par("usr")
		plot.width <- plot.dim[2] - plot.dim[1]
		plot.height <- plot.dim[4] - plot.dim[3]

		xmin <- plot.dim[2] - plot.width * 0.02
		xmax <- plot.dim[2] + plot.width * 0.02
		xmid <- mean(c(xmin, xmax))
		ymin <- plot.dim[3]
		ymax <- plot.dim[4]
		yseg <- segment.region(ymin, ymax, num.states+1, alignment = "ends")
		for(i in 1:(length(yseg)-1)){
			draw.rectangle(xmin, xmax, yseg[i], yseg[i+1], fill = state.col[i])
			text(x = xmid, y = mean(c(yseg[i], yseg[i+1])), labels = i)
		}

		strand <- unique(as.numeric(gene.info[,"strand"]))
	
		if(strand == -1){
			gene.start <- unique(as.numeric(gene.info[,"start_position"]))
			gene.end <- unique(as.numeric(gene.info[,"end_position"]))
			trans.dir <- "/\\"
		}else{
			gene.end <- unique(as.numeric(gene.info[,"start_position"]))
			gene.start <- unique(as.numeric(gene.info[,"end_position"]))
			# trans.dir <- "/\\"
			trans.dir <- "\\/"
		}

		gene.start.idx <- get.idx(rownames(state.mat), gene.start)
		gene.end.idx <- get.idx(rownames(state.mat), gene.end)
		
		fp.utr.start <- as.numeric(gene.info[,"5_utr_start"])
		fp.utr.start <- fp.utr.start[which(!is.na(fp.utr.start))]
		fp.utr.end <- as.numeric(gene.info[,"5_utr_end"])
		fp.utr.end <- fp.utr.end[which(!is.na(fp.utr.end))]
		
		tp.utr.start <- as.numeric(gene.info[,"3_utr_start"])
		tp.utr.start <- tp.utr.start[which(!is.na(tp.utr.start))]
		tp.utr.end <- as.numeric(gene.info[,"3_utr_end"])
		tp.utr.end <- tp.utr.end[which(!is.na(tp.utr.end))]
	
		tss <- sort(unique(as.numeric(gene.info[,"transcription_start_site"])))
		tss.idx <- unlist(lapply(tss, function(x) get.idx(rownames(state.mat), x)))
		
		if(is.finite(gene.end.idx) && is.finite(gene.start.idx)){
			arrows(y0 = 0.05, x1 = gene.start.idx, x0 = gene.end.idx, lwd = 2, length = 0.1)
		}else{
			if(strand == -1){
				text(y = 0.01, x = gene.start.idx, labels = trans.dir)
			}else{
				text(y = 0.01, x = gene.end.idx, labels = trans.dir)	
			}
		}

		
		#also mark the 5' and 3' UTRs, and TSS
		if(length(fp.utr.start) > 0){
			for(f in 1:length(fp.utr.start)){
				fps.idx <- get.idx(rownames(state.mat), fp.utr.start[f])
				fpe.idx <- get.idx(rownames(state.mat), fp.utr.end[f])
				if(is.finite(fps.idx) && is.finite(fpe.idx)){
					#segments(y0 = 0.14, x1 = fps.idx, x0 = fpe.idx, lwd = 3, col = "red")
				}
			}
		}
	
		if(length(tp.utr.start) > 0){
			for(f in 1:length(tp.utr.start)){
				tps.idx <- get.idx(rownames(state.mat), tp.utr.start[f])
				tpe.idx <- get.idx(rownames(state.mat), tp.utr.end[f])
				if(is.finite(tps.idx) && is.finite(tpe.idx)){
					#segments(y0 = 0.14, x1 = tps.idx, x0 = tpe.idx, lwd = 3, col = "blue")
				}
			}
		}


		par(xpd = FALSE)
		if(!is.null(gene.expr[[1]])){
			par(mar = c(4,4,2,2))
			stripchart(ordered.expr, vertical = TRUE, main = "Expression in Controls", 
			col = ordered.colors, pch = 16, cex = 1.5, axes = FALSE)
			axis(2)
			abline(h = median(unlist(ordered.expr)))
			plot.dim <- par("usr")
			plot.height <- plot.dim[4]-plot.dim[3]
			box.labels <- unlist(lapply(strsplit(names(ordered.expr), "/"), function(x) x[1]))
			par(xpd = TRUE)
			text(x = 1:length(ordered.expr), y = (plot.dim[3]-(plot.height*0.1)), 
			labels = box.labels, adj = 0.5)
			par(xpd = FALSE)

			#if state.weights are supplied, use these to calculate the
			#correlation with expression. Otherwise, use the specified
			#dimension reduction strategy
			if(!is.null(state.weights)){
				u_states <- 1:num.states
				prop.mat <- matrix(NA, ncol = ncol(ordered.states), nrow = length(u_states))
				colnames(prop.mat) <- colnames(ordered.states)
				rownames(prop.mat) <- u_states
				for(st in u_states){
					prop.mat[st,] <- apply(ordered.states, 2, 
					function(x) length(which(x == st))/length(x))
				}
				weighted.prop <- apply(prop.mat, 2, function(x) x*state.weights)
				scaled.mat <- matrix(colSums(weighted.prop), ncol = 1)
				ytext <- "Weighted States"
			}else{
				if(dim.rd == "mds"){
					state.dist <- get.hamming.mat(ordered.states)
					scaled.mat <- cmdscale(state.dist, k = 1)
					ytext <- "Scaled State"
				}else{
					scaled.mat <- matrix(apply(ordered.states, 2, 
					function(x) length(which(x == state))/length(x)), ncol = 1)
					ytext <- paste0("Proportion State", state)
				}
			}	

			par(mar = c(7,4,2,2))
			barplot(scaled.mat[,1], col = ordered.colors, ylab = ytext,  
			names.arg = ordered.names, las = 2)
			mean.expr <- sapply(ordered.expr, mean)

			scaled.expr <- scale(mean.expr)
			model <- lm(scaled.expr~scaled.mat[,1])
			expr.cor <- cor.test(mean.expr, scaled.mat[,1], use = "complete")
			r <- signif(expr.cor$estimate, 2)
			p <- signif(expr.cor$p.value, 2)
			ext.name <- gene.info$external_gene_name[1]
			
			plot(scaled.mat[,1], scaled.expr, col = ordered.colors, pch = 16, 
			xlab = ytext, ylab = paste(ext.name, "expression"), cex = 1.5, 
			main = paste0(ext.name, ": r = ", r, ", p = ", p))
			abline(model)
		}
	
	mtext(paste(unique(gene.info[,"external_gene_name"]), " Chr", 
	unique(gene.info[,"chromosome_name"]), ":", unique(gene.info[,"start_position"]), "..", 
	unique(gene.info[,"end_position"]), sep = ""), outer = TRUE, line = -2)
	
	}
