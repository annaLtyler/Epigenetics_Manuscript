#This gene gets the ChromHMM state proportions of 
#a given gene
get.chrom.state <- function(bed.info, rnaseq.gene.info, rna.seq, gene.id, start.feature, end.feature, upstream.buffer, downstream.buffer){
	
	num.strains <- length(bed.info)
	strain.names <- names(bed.info)
	gene.locale <- which(rnaseq.gene.info[,1] == gene.id)
	gene.expr.locale <- which(rownames(rna.seq) == gene.id)
	
	if(length(gene.locale) > 0){
		gene.info <- rnaseq.gene.info[gene.locale,,drop=FALSE]
	
		#use larger buffers then trim to make sure we get states
		strain.states <- lapply(bed.info, function(x) get.states(gene.info, x, upstream.buffer = 50000, downstream.buffer = 50000, start.feature = start.feature, end.feature = end.feature))
		names(strain.states) <- strain.names
			
		if(sum(unlist(lapply(strain.states, nrow))) > 0){
			
			#line up the states for the different strains and trim to requested region
			gene.mat <- align.states(strain.states)
			gene.mat <- trim.state.mat(state.mat = gene.mat, gene.info, upstream.buffer, 
			downstream.buffer, start.feature, end.feature)
			}else{
			gene.mat <- NA	
			}	
		return(gene.mat)
		}else{
		return(NA)
		}

	
	}

