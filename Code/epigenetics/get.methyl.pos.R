get.methyl.pos <- function(all.methyl.data, strain.name = NULL, treatment.name, 
id.table, chr, start.pos, stop.pos, average.replicates = TRUE){

	if(is.null(strain.name)){strain.name <- unique(id.table[,"strain"])}

	strain.locale <- which(id.table[,1] %in% strain.name)
	treatment.locale <- which(id.table[,2] %in% treatment.name)
	id.locale <- intersect(strain.locale, treatment.locale)
	
	methyl.list <- vector(mode = "list", length = length(id.locale))
	names(methyl.list) <- apply(id.table[id.locale, ], 1, function(x) paste(x, collapse = "_"))
	
	
	for(i in 1:length(id.locale)){
		methyl.table <- all.methyl.data[[id.locale[i]]]
		chr.locale <- which(methyl.table[,1] == chr)
		chr.table <- methyl.table[chr.locale,]
		after.start <- which(chr.table[,2] >= start.pos)
		before.end <- which(chr.table[,3] <= stop.pos)
		section.locale <- intersect(after.start, before.end)
		if(length(section.locale) == 0){return("No methylation data in this range.")}
		if(length(section.locale) > 0){
			methyl.list[[i]] <- chr.table[section.locale,2:4]
			}
		}
	
	has.entries <- length(which(unlist(lapply(methyl.list, length)) > 0))
	
	if(has.entries > 0){
		methyl.mat <- align.states(state.list = methyl.list, chunk.size = 1, remove.nas = TRUE)
		non.rep.ids <- colnames(methyl.mat)
		if(average.replicates){
			non.rep.ids <- apply(id.table[id.locale,], 1, function(x) paste(x[1:2], collapse = "_"))
			#average over replicates
			u_ids <- unique(non.rep.ids)
			strain.table <- matrix(NA, nrow = nrow(methyl.mat), ncol = length(u_ids))
			for(i in 1:length(u_ids)){
				# strain.locale <- which(non.rep.ids[id.locale] == u_ids[i])
				strain.locale <- which(non.rep.ids== u_ids[i])
				strain.table[,i] <- rowMeans(methyl.mat[,strain.locale,drop=FALSE], na.rm = TRUE)
				}
			}else{
			strain.table <- methyl.mat
			u_ids <- non.rep.ids
			}
		colnames(strain.table) <- u_ids
		rownames(strain.table) <- rownames(methyl.mat)
		
		return(strain.table)
		# imageWithText(t(strain.table), show.text = FALSE, row.names = u_ids, split.at.vals = TRUE, split.points = 50)
		}
		
	
}