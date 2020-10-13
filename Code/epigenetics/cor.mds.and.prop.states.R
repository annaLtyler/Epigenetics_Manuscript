cor.mds.and.prop.states <- function(chrom.mats, chrom.prop.mats, name.color.key){
	
	
	chrom.strain.order <- match(colnames(chrom.mats[[1]]), name.color.key[,4])
	prop.strain.order <- match(colnames(chrom.prop.mats[[1]]), name.color.key[,4])
	
	chrom.strains <- name.color.key[chrom.strain.order,1]
	prop.strains <- name.color.key[prop.strain.order,1]
	shared.strains <- intersect(chrom.strains, expr.strains)
	
	prop.strain.locale <- match(shared.strains, prop.strains)
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
	
	all.cor <- matrix(NA, nrow = length(scaled.v), ncol = nrow(chrom.prop.mats[[1]]))
	rownames(all.cor) <- names(scaled.v)
	colnames(all.cor) <- 1:ncol(all.cor)
	for(i in 1:length(scaled.v)){
		report.progress(i, length(scaled.v))
		gene.id <- names(scaled.v)[i]
		chrom.prop.locale <- which(names(chrom.prop.mats) == gene.id)
		all.cor[i,] <- apply(chrom.prop.mats[[chrom.prop.locale]], 1, function(x) cor(x[prop.strain.locale], scaled.v[[i]][chrom.strain.locale]))
		}
	return(all.cor)	
}