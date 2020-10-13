#This function plots eqtl results
	
plot.eqtl <- function(coef.list, lod.score, label, gene.coord = NULL, max.lod = NULL, min.coef = NULL, max.coef = NULL){
		
	do.cols <- get.allele.colors()
	
	layout.mat <- matrix(c(1,2,0), ncol = 1, byrow = TRUE)
	layout(layout.mat, heights = c(0.5, 1, 0.3))
	
	coords <- as.numeric(unlist(lapply(strsplit(rownames(coef.list), "_"), function(x) x[2])))
	
	par(mar = c(0,4,4,4))
	
	if(is.null(max.lod)){
		max.lod <- max(lod.score, na.rm = TRUE)
		}
	plot.height.lod <- max(lod.score, na.rm = TRUE) - min(lod.score, na.rm = TRUE)
	plot(x = coords, y = lod.score, ylim = c(0, max.lod), type = "l", axes = FALSE, xlab = "", ylab = "LOD score")
	abline(v = gene.coord)
	# points(x = gene.coord, y = 0, pch = "*", col = "red", cex = 2)
	axis(2)	
	par(xpd = TRUE)
	text(x = median(coords), y = max.lod+(plot.height.lod*0.15), labels = label, cex = 1.3)
	par(xpd = FALSE)
	
	par(mar = c(6,4,0,4))
	if(is.null(max.coef)){
		max.coef <- (max(coef.list, na.rm = TRUE)+(max(coef.list, na.rm = TRUE)*0.1))
		}
	if(is.null(min.coef)){
		min.coef <- min(coef.list, na.rm = TRUE)
		}
		plot.height <- max.coef - min.coef
	plot.new()
	plot.window(xlim = c(min(coords), max(coords)), ylim = c(min.coef, max.coef))

	for(i in 1:ncol(coef.list)){
		strain.locale <- which(do.cols[,2] == colnames(coef.list)[i])
		points(x = coords, y = coef.list[,i], col = do.cols[strain.locale,3], type = "l", lwd = 2)
		}
	abline(v = gene.coord)
	axis(1);axis(2)
	par(xpd = TRUE)
	legend(x = 1, y = min.coef-(plot.height*0.2), legend = colnames(coef.list), col = do.cols[,3], lty = 1, lwd = 2, horiz = TRUE, cex = 0.7)
	par(xpd = FALSE)

}

