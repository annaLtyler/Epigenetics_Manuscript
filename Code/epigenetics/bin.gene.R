#This function creates biologically informed gene bins
#upstream, first exon, other exons, introns, and downstream

bin.gene <- function(states, strand, tss, gene.end, exon.table, collapse.exons = TRUE){

	bin.buffer <- 0
	
	if(strand == -1){
		new.start <- gene.end
		new.end <- tss
		tss <- new.start
		gene.end <- new.end
		}
		
	if(strand == 1){
	upstream.bins <- which(rownames(states) <= tss)
	downstream.bins <- which(rownames(states) >= gene.end)
	promoter.bins <- intersect(which(rownames(states) >= tss - (700+bin.buffer)), which(rownames(states) <= tss + (200+bin.buffer)))
	}else{
	upstream.bins <- which(rownames(states) >= tss)
	downstream.bins <- which(rownames(states) <= gene.end)
	promoter.bins <- intersect(which(rownames(states) <= tss + (700+bin.buffer)), which(rownames(states) >= tss - (200+bin.buffer)))
	}
	
	
	get.exon.bins <- function(exon.start, exon.end){
		after.start <- which(rownames(states) >= (as.numeric(exon.start)-bin.buffer))
		before.end <- which(rownames(states) <= (as.numeric(exon.end)+bin.buffer))		
		return(intersect(after.start, before.end))
		}

	get.intron.bins <- function(exon.end, exon.start){
		after.end <- which(rownames(states) >= (as.numeric(exon.end)-bin.buffer))
		before.start <- which(rownames(states) <= (as.numeric(exon.start)+bin.buffer))		
		return(intersect(after.end, before.start))
		}


	all.exon.bins <- apply(exon.table, 1, function(x) get.exon.bins(x[3], x[4]))
	if(collapse.exons){
		all.exon.bins <- as.vector(unlist(all.exon.bins))
		}else{
		names(all.exon.bins) <- paste0("exon", 1:length(all.exon.bins))
		}

	all.intron.bins <- vector(mode = "list", length = nrow(exon.table))

	if(nrow(exon.table) > 1){
		for(i in 1:(nrow(exon.table)-1)){
			all.intron.bins[[i]] <- get.intron.bins(exon.table[i,4], exon.table[(i+1), 3])
			}
		}
	if(collapse.exons){
		all.intron.bins <- as.vector(unlist(all.intron.bins))
		}else{
		names(all.intron.bins) <- paste0("intron", 1:length(all.intron.bins))
		}
	
	if(collapse.exons){
		all.bins <- list("upstream" = upstream.bins, "promoter" = promoter.bins, "downstream" = downstream.bins, "exon" = all.exon.bins, "intron" = all.intron.bins)
	}else{
	all.bins <- list("upstream" = upstream.bins, "promoter" = promoter.bins, "downstream" = downstream.bins)
	all.bins <- c(all.bins, all.exon.bins, all.intron.bins)
	}
	return(all.bins)

}