#This gene gets the ChromHMM state proportions of 
#up and downstream of a given point

get.chrom.state.nearby <- function(bed.info, chr, loc, upstream.buffer, downstream.buffer, num.states = 6){
	
	
	num.strains <- length(bed.info)
	strain.names <- names(bed.info)
			
	chr.locale <- lapply(bed.info, function(x) which(x[,1] == paste0("chr", chr)))		
	after.start <- lapply(bed.info, function(x) which(x[,2] >= loc-upstream.buffer))
	before.end <- lapply(bed.info, function(x) which(x[,3] <= loc+downstream.buffer))
	
	in.gene <- vector(mode = "list", length = length(chr.locale))
	names(in.gene) <- strain.names
	for(i in 1:length(chr.locale)){
		in.gene[[i]] <- Reduce(intersect, list(chr.locale[[i]], after.start[[i]], before.end[[i]]))			
		}		
	
	gene.states <- vector(mode = "list", length = length(chr.locale))
	names(gene.states) <- strain.names
	for(i in 1:length(gene.states)){
		gene.states[[i]] <- bed.info[[i]][in.gene[[i]],2:4]	
		}
				
		if(sum(unlist(lapply(gene.states, length))) > 0){
			
			#line up the states for the different strains and trim to requested region
			gene.mat <- align.states(gene.states)

			#remove the rows with NAs
			na.idx <- which(is.na(gene.mat), arr.ind = TRUE)
			if(nrow(na.idx) > 0){
				gene.mat <- gene.mat[-unique(na.idx[,1]),]
				}
			}else{
			gene.mat <- NA	
			}	
		return(gene.mat)
		

	}

