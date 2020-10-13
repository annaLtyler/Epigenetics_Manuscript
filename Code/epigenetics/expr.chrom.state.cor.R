#correlate mean inbred expression to scaled chromatin
#matrix

expr.chrom.state.cor <- function(chrom.mats, group.mean.expr, name.color.key){

	
	chrom.strain.order <- match(colnames(chrom.mats[[1]]), name.color.key[,4])
	expr.strain.order <- match(names(group.mean.expr[[1]]), name.color.key[,4])
	chrom.strains <- name.color.key[chrom.strain.order,1]
	expr.strains <- name.color.key[expr.strain.order,1]
	shared.strains <- intersect(chrom.strains, expr.strains)
	
	expr.strain.locale <- match(shared.strains, expr.strains)
	chrom.strain.locale <- match(shared.strains, chrom.strains)
	
	is.good <- function(mat){
		if(is.null(mat)){
			return(FALSE)
			}
		if(is.na(mat)[1]){
			return(FALSE)
			}
		if(is.null(dim(mat))){
			return(FALSE)
			}
		return(TRUE)
		}
	good.locale <- which(unlist(lapply(chrom.mats, is.good)))
	
	
	good.mats <- chrom.mats[good.locale]
	cat("Calculating Hamming Distance...\n")
	all.hamming <- lapply(good.mats, function(x) get.hamming.mat(x))	
	names(all.hamming) <- names(chrom.mats)[good.locale]
	
	good.hamm <- which(unlist(lapply(all.hamming, length)) > 0)
	
	cat("Scaling Matrices...\n")
	scaled.v <- lapply(all.hamming[good.hamm], function(x) cmdscale(x, 1))
	
	all.cor <- rep(NA, length(good.hamm))
	all.p <- rep(NA, length(good.hamm))
	
	for(i in 1:length(good.hamm)){
		report.progress(i, length(good.hamm))
		gene.id <- names(scaled.v)[good.hamm[i]]
		expr.locale <- which(names(group.mean.expr) == gene.id)
		if(length(expr.locale) > 0){
			if(ncol(scaled.v[[i]]) > 0 && !is.na(group.mean.expr[[expr.locale]][1])){
				# plot(scaled.v[[i]][chrom.strain.locale,], group.mean.expr[[expr.locale]][expr.strain.locale])
				result <- cor.test(scaled.v[[i]][chrom.strain.locale,], group.mean.expr[[expr.locale]][expr.strain.locale])
				all.cor[i] <- result$estimate
				all.p[i] <- result$p.value
				}
			}
		}
	
	final.table <- cbind(all.cor, all.p)
	colnames(final.table) <- c("Correlation", "P.value")
	rownames(final.table) <- names(all.hamming)[good.hamm]
	return(final.table)
	}