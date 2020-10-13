#This function looks at the correlation between 
#the chromatin state proportions multiplied by the 
#state code and expression, as well as the correlation
#between eQTL coefficients and expression. 
# col = "#de2d26" #red
# col = "#3182bd" #blue

Chromatin.Coef.Correlation <- function(do.eqtl, coef.list, chrom.states, state.code, 
group.mean.expr, lod.list = c(10, 20, 50), cut.point = 0.5, col = "#de2d26", 
use.state.code = TRUE, num.strains = 9, start.feature, end.feature, upstream.buffer, 
downstream.buffer, path = "."){

if(use.state.code){
	jpg.text <- "state.code"
	}else{
	jpg.text <- "scaled.matrix"	
	}
pear.chrom.list <- vector(mode = "list", length = length(lod.list))
names(pear.chrom.list) <- paste("LOD_", lod.list, sep = "")

max.lod <- unlist(lapply(do.eqtl, function(x) max(x[[3]])))
do.eqtl.gene.name <- unlist(lapply(do.eqtl, function(x) x[[1]][1,2]))
do.eqtl.gene.id <- unlist(lapply(do.eqtl, function(x) x[[1]][1,1]))


eqtl.info <- function(do.eqtl.ind){
	max.score <- max(do.eqtl.ind[[3]])
	gene.id <- do.eqtl.ind[[1]][1,1]
	gene.name <- do.eqtl.ind[[1]][1,2]
	return(c(gene.name, gene.id, max.score))
	}

gene.info <- lapply(do.eqtl, eqtl.info)

for(l in 1:length(lod.list)){

	filename <- paste("Correlation.Between.Inbred.Expr.and.DO.coef.LOD.", lod.list[l], ".", 
	jpg.text, ".", start.feature, ".", end.feature, ".", upstream.buffer, ".", downstream.buffer, 
	".jpg", sep = "")
	jpeg(file.path(path, filename), width = 10, height = 10, units = "in", res = 300)
	par(mfrow = c(2,2))

	cat("LOD:", lod.list[l], "\n")
	good.lod.locale <- which(max.lod >= lod.list[l])
	all.pear.cor <- rep(NA, length(good.lod.locale))
	names(all.pear.cor) <- do.eqtl.gene.name[good.lod.locale]
	all.pear.cor.chrom <- rep(NA, length(good.lod.locale))
	names(all.pear.cor.chrom) <- do.eqtl.gene.name[good.lod.locale]
	for(i in 1:length(good.lod.locale)){
		report.progress(i, length(good.lod.locale))
		gene.id <- do.eqtl[good.lod.locale[i]][[1]][[1]][1,1]
		gene.locale <- which(names(coef.list) == gene.id)
		cis.coef <- coef.list[[gene.locale]]
		if(var(cis.coef) != 0){
			coef.order <- match(col.table[,5], names(cis.coef))

			#look at correlation with expression
			gene.locale <- which(names(group.mean.expr) == gene.id)
			if(length(gene.locale) > 0){
				inbred.expr <- group.mean.expr[[gene.locale]]
				if(!is.na(inbred.expr[1])){
					expr.order <- match(col.table[,4], names(inbred.expr))
					all.pear.cor[i] <- cor(cis.coef[coef.order], 
					inbred.expr[expr.order], use = "complete.obs", method = "pearson")
					# plot(cis.coef[coef.order], inbred.expr[expr.order])
					}
				}
				
				#and correlation with chromatin
				if(use.state.code){
				gene.locale <- which(names(chrom.states) == gene.id)
				if(length(gene.locale) > 0){
					gene.states <- chrom.states[[gene.locale]]
					if(!is.na(gene.states) && sum(gene.states) > 0){
						coded.states <- apply(gene.states, 2, function(x) x*state.code)
						weighted.states <- colMeans(coded.states, na.rm = TRUE)
						all.pear.cor.chrom[i] <- cor(weighted.states, 
						inbred.expr, method = "pearson")
						#plot(weighted.states, inbred.expr)
						}
					}
				}else{#end use.state.code
				gene.locale <- which(names(chrom.mats) == gene.id)
				if(length(gene.locale) > 0){
					gene.mat <- chrom.mats[[gene.locale]]
					if(class(gene.mat) != "matrix"){
						gene.mat <- matrix(gene.mat, ncol = num.strains)
						colnames(gene.mat) <- colnames(chrom.mats[[gene.locale]])
						}
					if(!is.na(gene.mat[1])){
						mat.order <- match(col.table[,6], colnames(gene.mat))
						d.mat <- get.hamming.mat(gene.mat[,mat.order])
						if(sum(d.mat) > 0){
							scaled.vector <- cmdscale(d.mat, 1)
							all.pear.cor.chrom[i] <- cor(inbred.expr[expr.order], scaled.vector[,1])
							#plot(inbred.expr[expr.order], scaled.vector[,1])
							}
						}
					}
				}#end use full chromatin matrix

			} 
		}#end looping through good locale
	
	hist(all.pear.cor, breaks = 100, main = paste("Expression\nMinimum LOD:", lod.list[l]), 
	xlab = "Pearson Correlation", xlim = c(-1, 1))
	if(use.state.code){
		hist(all.pear.cor.chrom, breaks = 100, main = paste("Chromatin\nMinimum LOD:", 
		lod.list[l]), xlab = "Pearson Correlation", xlim = c(-1, 1))
		}else{
		hist(abs(all.pear.cor.chrom), breaks = 100, main = paste("Chromatin\nMinimum LOD:", 
		lod.list[l]), xlab = "Pearson Correlation", xlim = c(0, 1))
			}
	
	if(use.state.code){
		plot(all.pear.cor, all.pear.cor.chrom, xlab = "eQTL correlation", 
		ylab = "Chromatin correlation")
		}else{
		plot(all.pear.cor, abs(all.pear.cor.chrom), xlab = "eQTL correlation", 
		ylab = "Chromatin correlation")	
		}
	abline(0,1, col = "red", lwd = 2)
	
	
	#calculate eQTL-Chromatin correlations for eQTL-expression correlations > 0.5
	high.cor <- which(all.pear.cor >= cut.point)
	low.cor <- which(all.pear.cor < cut.point)
	high.chrom.cor <- which(all.pear.cor.chrom >= cut.point)
	low.chrom.cor <- which(all.pear.cor.chrom < cut.point)
	high.high <- length(intersect(high.cor, high.chrom.cor))
	
	
	high.high.genes <- all.pear.cor[intersect(high.cor, high.chrom.cor)]
	high.high.genes.chrom <- all.pear.cor.chrom[intersect(high.cor, high.chrom.cor)]
	high.high.table <- cbind(high.high.genes, high.high.genes.chrom)
	high.high.table <- sort.by.then.by(high.high.table, c(1,2), c("n", "n"), decreasing = TRUE)
	# plot.gene.list.stats(rev(rownames(high.high.table)), "high.chromatin.high.eQTL.cor")
	
	high.low <- length(intersect(high.cor, low.chrom.cor))
	low.high <- length(intersect(low.cor, high.chrom.cor))
	low.low <- length(intersect(low.cor, low.chrom.cor))
	
	if(use.state.code){
	plot(all.pear.cor, all.pear.cor.chrom, xlab = "eQTL correlation", 
	ylab = "Chromatin correlation", col = rgb(0,0,0, alpha = 0.2))
	cor.min <- -1
	}else{
	plot(all.pear.cor, abs(all.pear.cor.chrom), xlab = "eQTL correlation", 
	ylab = "Chromatin correlation", col = rgb(0,0,0, alpha = 0.2))		
	cor.min = 0
		}
	abline(h = cut.point, col = col, lwd = 2);abline(v = cut.point, col = col, lwd = 2)
	text(x = mean(c(cut.point,1)), y = mean(c(cut.point,1)), labels = high.high, 
	col = col, cex = 2, font = 2)
	text(x = mean(c(cut.point,-1)), y = mean(c(cut.point,cor.min)), labels = low.low, 
	col = col, cex = 2, font = 2)
	text(x = mean(c(cut.point,1)), y = mean(c(cut.point,cor.min)), labels = high.low, 
	col = col, cex = 2, font = 2)
	text(x = mean(c(cut.point,-1)), y = mean(c(cut.point,1)), labels = low.high, 
	col = col, cex = 2, font = 2)
	
	#calculate number of correlations at each 0.5 percentage level
	percentage.seq <- seq(-1, 1, 0.05)
	perc.above <- rep(NA, length(percentage.seq))
	perc.above.chrom <- rep(NA, length(percentage.seq))
	num.above <- rep(NA, length(percentage.seq))
	num.above.chrom <- rep(NA, length(percentage.seq))
	
	for(p in 1:length(percentage.seq)){
		perc.above[p] <- length(which(all.pear.cor >= percentage.seq[p]))/length(which(!is.na(all.pear.cor)))
		perc.above.chrom[p] <- length(which(all.pear.cor.chrom >= percentage.seq[p]))/length(which(!is.na(all.pear.cor.chrom)))
		num.above[p] <- length(which(all.pear.cor >= percentage.seq[p]))
		num.above.chrom[p] <- length(which(all.pear.cor.chrom >= percentage.seq[p]))
		}
	final.table <- cbind(signif(percentage.seq, 2), num.above, round(perc.above*100), 
	num.above.chrom, round(perc.above.chrom*100))
	colnames(final.table) <- c("Pearson.Correlation", "Num.Above.eQTL", "Percent.Above.eQTL", 
	"Num.Above.Chrom", "Percent.Above.Chrom")
	
	filename <- paste("Percentage.Above.Correlation.LOD.", lod.list[l], ".", 
	start.feature, ".", end.feature, ".", upstream.buffer, ".", downstream.buffer, ".txt", 
	sep = "")

	cat("\nWriting to", filename, '\n')
	
	write.table(final.table, file.path(path, filename), quote = FALSE, sep = "\t", 
	row.names = FALSE)

	result.table <- cbind(all.pear.cor, all.pear.cor.chrom)
	colnames(result.table) <- c("Coef-Expr", "Chrom-Expr")
	na.locale <- unique(which(is.na(result.table), arr.ind = TRUE)[,1])
	result.table <- result.table[-na.locale,]
	
	filename <- paste("Gene.Cor.Table.LOD.", lod.list[l], ".", start.feature, ".", 
	end.feature, ".", upstream.buffer, ".", downstream.buffer, ".txt", sep = "")
	
	cat("\nWriting to", filename, "\n")
	write.table(result.table, file.path(path, filename), quote = FALSE, sep = "\t")
	pear.chrom.list[[l]] <- result.table
	dev.off()
} #end looping through LOD scores


	return(pear.chrom.list)
}