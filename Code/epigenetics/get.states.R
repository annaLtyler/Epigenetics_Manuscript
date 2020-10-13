#This function retrieves ChromHMM states around a given gene
#state info is the bed file containing the ChromHMM states
#gene.info is the output of get.gene.info. It contains position
#and name information about a gene
get.states <- function(gene.info, state.info, upstream.buffer, downstream.buffer, start.feature = "start_position", end.feature = "end_position"){

	strand <- unique(gene.info[,"strand"])
	
	if(strand == -1){
		feature.names <- switch.strand(start.feature, end.feature)
		start.feature <- feature.names[1]
		end.feature <- feature.names[2]
		}
		
	gene.chr <- unique(gene.info[,"chromosome_name"])
	chr.locale <- which(state.info[,1] == paste("chr", gene.chr, sep = ""))

	gene.start <- unique(gene.info[,start.feature])
	gene.start <- gene.start[which(!is.na(gene.start))]
	gene.end <- unique(gene.info[,end.feature])
	gene.end <- gene.end[which(!is.na(gene.end))]
	
	if(strand == 1){
		after.start <- which(as.numeric(state.info[,2]) >= as.numeric(gene.start - upstream.buffer))
		before.end <- which(as.numeric(state.info[,3]) <= as.numeric(gene.end + downstream.buffer))
		}else{
		after.start <- which(as.numeric(state.info[,2]) <= as.numeric(gene.start + upstream.buffer))
		before.end <- which(as.numeric(state.info[,3]) >= as.numeric(gene.end - downstream.buffer))
		}
	
	in.gene <- Reduce(intersect, list(chr.locale, after.start, before.end))
	gene.states <- state.info[in.gene,c(2:4)]
	return(gene.states)
	}
