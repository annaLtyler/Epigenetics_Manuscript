#Given a list of states with their positions, this 
#function aligns the states into a single matrix
#so we can compare states across strains for example

align.methylation <- function(methyl.list, strain.means = TRUE){
	require(stringr)
	
	chunk.size = 1
	chr <- unique(unlist(lapply(methyl.list, function(x) x[,1])))
	chr.locale <- lapply(methyl.list, function(x) which(x[,1] == chr))
		
	min.start <- rep(NA, length(methyl.list))
	max.end <- rep(NA, length(methyl.list))
		
		for(m in 1:length(methyl.list)){
			min.start[m] <- min(methyl.list[[m]][chr.locale[[m]],2])
			max.end[m] <- max(methyl.list[[m]][chr.locale[[m]],3])
			}
	
	gene.seq <- seq(min(min.start), max(max.end), chunk.size)
	state.mat <- matrix(NA, nrow = length(gene.seq), ncol = length(methyl.list))

	
		for(i in 1:length(methyl.list)){
			start.locale <- match(methyl.list[[i]][,2], gene.seq)
			end.locale <- match(methyl.list[[i]][,3], gene.seq)
			if(length(start.locale) > 0){
				for(j in 1:length(start.locale)){
					state.mat[start.locale[j]:end.locale[j],i] <- methyl.list[[i]][j,4]
					}
				}
			}
			
		
		rownames(state.mat) <- gene.seq
		colnames(state.mat) <- unlist(lapply(strsplit(names(methyl.list), "_"), function(x) x[3]))

		all.na <- apply(state.mat, 1, function(x) all(is.na(x)))
		state.mat <- state.mat[which(!all.na),]

		# imageWithText(t(state.mat), show.text = FALSE, row.names = colnames(state.mat))
		colnames(state.mat) <- str_to_upper(colnames(state.mat))
	
		if(strain.means){
			u_strains <- unique(colnames(state.mat))
			means <- lapply(u_strains, function(x) rowMeans(state.mat[,which(colnames(state.mat) == x),drop=FALSE], na.rm = TRUE))
			means.mat <- Reduce("cbind", means)
			colnames(means.mat) <- u_strains
			state.mat <- means.mat
			}

		# imageWithText(t(state.mat), show.text = FALSE, row.names = colnames(state.mat))

		return(state.mat)
	}
