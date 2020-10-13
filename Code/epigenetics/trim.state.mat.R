#This function trims a state matrix to be around 
#given gene features +/- buffers

trim.state.mat <- function(state.mat, gene.info, upstream.buffer, downstream.buffer, start.feature, end.feature){

	strand <- unique(gene.info[,"strand"])
	
	if(strand == -1){
		feature.names <- switch.strand(start.feature, end.feature)
		start.feature <- feature.names[1]
		end.feature <- feature.names[2]
		}

	gene.start <- unique(gene.info[,start.feature])
	gene.start <- min(gene.start[which(!is.na(gene.start))])
	gene.end <- unique(gene.info[,end.feature])
	gene.end <- max(gene.end[which(!is.na(gene.end))])

	if(strand == 1){
		start.idx <- get.idx(as.numeric(rownames(state.mat)), (as.numeric(gene.start) - upstream.buffer))
		# rownames(state.mat)[start.idx]
		end.idx <- get.idx(as.numeric(rownames(state.mat)), (as.numeric(gene.end) + downstream.buffer))
		# rownames(state.mat)[end.idx]
		if(!is.finite(start.idx)){start.idx <- 1}
		if(!is.finite(end.idx)){end.idx <- nrow(state.mat)}
		good.vals <- c(start.idx:end.idx)
		# pos.vector <- c("gene.start" = gene.start, "gene.end" = gene.end, "state.start" = as.numeric(rownames(state.mat)[start.idx]), "state.end" = as.numeric(rownames(state.mat)[end.idx]))
		# sorted.vector <- sort(pos.vector)
		# sorted.vector
		}else{
		start.idx <- get.idx(as.numeric(rownames(state.mat)), (as.numeric(gene.start) + upstream.buffer))
		end.idx <- get.idx(as.numeric(rownames(state.mat)), (as.numeric(gene.end) - upstream.buffer))
		if(!is.finite(start.idx)){start.idx <- nrow(state.mat)}
		if(!is.finite(end.idx)){end.idx <- 1}
		good.vals <- c(end.idx:start.idx)
		}
	
	
	new.state.mat <- state.mat[good.vals,]
	return(new.state.mat)
	
	}
