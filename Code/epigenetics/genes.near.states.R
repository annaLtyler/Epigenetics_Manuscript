genes.near.states <- function(all.bed, state.num){	
	
	chr = unique(all.bed[[1]][,1])
	chr <- chr[which(chr != "chrY")]
	all.dist.mats <- vector(mode = "list", length = length(all.bed))
	names(all.dist.mats) <- names(all.bed)

	for(i in 1:length(chr)){
		cat(chr[i], "\n")
		
		chr.locale <- lapply(all.bed, function(x) which(x[,1] == chr[i]))
		state.locale <- lapply(all.bed, function(x) which(x[,4] == state.num))
		chr.state.overlap <- vector(mode = "list", length = length(chr.locale))
		for(j in 1:length(all.bed)){
			chr.state.overlap[[j]] <- intersect(chr.locale[[j]], state.locale[[j]])
			}

		for(j in 1:length(all.bed)){
			cat("\t", names(all.bed)[j], "\n")
			state.dist.to.nearest.gene <- lapply(chr.state.overlap[[j]], function(x) dist.to.nearest.TSS(gene.info.table, chr = gsub("chr", "", chr[i]), start.pos = all.bed[[j]][x,2], end.pos = all.bed[[j]][x,3]))
			
			all.dist.mats[[j]] <- rbind(all.dist.mats[[j]], Reduce(rbind, state.dist.to.nearest.gene))
			}
	
	} #end looping through Chr
	return(all.dist.mats)
	
}