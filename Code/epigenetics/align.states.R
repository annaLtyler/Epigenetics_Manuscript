#Given a list of states with their positions, this 
#function aligns the states into a single matrix
#so we can compare states across strains for example

align.states <- function(state.list, chunk.size = 200, remove.nas = FALSE){

	min.start <- min(unlist(lapply(state.list, function(x) if(nrow(x) > 0 && length(x) > 0){min(x[,1])})))
	max.end <- max(unlist(lapply(state.list, function(x) if(nrow(x) > 0 && length(x) > 0){max(x[,2])})))
	gene.seq <- seq(min.start, max.end, chunk.size)
	state.mat <- matrix(NA, nrow = length(gene.seq), ncol = length(state.list))
	
	for(i in 1:length(state.list)){
		#find the start and stop of each state in the state.list
		start.locale <- match(state.list[[i]][,1], gene.seq)
		end.locale <- match(state.list[[i]][,2], gene.seq)
		if(length(start.locale) > 0){
			for(j in 1:length(start.locale)){
				if(end.locale[j] - start.locale[j] == 1){
					state.mat[start.locale[j],i] <- state.list[[i]][j,3]
					}else{
					state.mat[start.locale[j]:end.locale[j],i] <- state.list[[i]][j,3]		
					}			
				}
			}
		}
	rownames(state.mat) <- gene.seq
	colnames(state.mat) <- names(state.list)
	if(remove.nas){
		all.na <- apply(state.mat, 1, function(x) length(which(is.na(x))))
		not.all.na <- which(all.na < ncol(state.mat))
		state.mat <- state.mat[not.all.na,]
		}
	# image(state.mat)
	return(state.mat)
}
