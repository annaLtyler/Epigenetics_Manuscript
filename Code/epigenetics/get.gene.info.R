#This function retrieves gene information from rnaseq.gene.info
#from epigen.states.and.expr

get.gene.info <- function(rnaseq.gene.info, rna.seq, all.bed, 
gene.id, T.or.C, start.feature, end.feature, upstream.buffer, downstream.buffer, 
get.snps = FALSE){
	
	gene.locale <- which(rnaseq.gene.info[,1] == gene.id)
	gene.expr.locale <- which(rownames(rna.seq) == gene.id)
	
	if(length(gene.locale) > 0){
		gene.info <- rnaseq.gene.info[gene.locale,,drop=FALSE]
	
		#use larger buffers then trim to make sure we get states
		strain.states <- lapply(all.bed, function(x) 
		get.states(gene.info, x, upstream.buffer = 50000, downstream.buffer = 50000, 
		start.feature = start.feature, end.feature = end.feature))
			
		if(sum(unlist(lapply(strain.states, nrow))) > 0){

			if(length(gene.expr.locale) > 0){
				#pull out expression values for controls
				gene.expr <- rna.seq[gene.expr.locale,]	
				labels <- substr(names(gene.expr), 1, 3)
				u_labels <- sort(unique(labels))
				binned.expr <- lapply(u_labels, function(x) as.vector(gene.expr[which(labels == x)]))
				names(binned.expr) <- u_labels
				treat.locale <- which(substr(names(binned.expr), 3, 3) == T.or.C)
				treat.expr <- binned.expr[treat.locale]
				names(treat.expr) <- substr(names(treat.expr), 1, 2)
				med.order <- order(unlist(lapply(binned.expr[treat.locale], median)))
			}else{
				med.order <- 1:ncol(rna.seq)	
			}

			
			#line up the states for the different strains and trim to requested region
			gene.mat <- align.states(strain.states)
			#imageWithText(t(gene.mat), show.text = FALSE, use.pheatmap.colors = TRUE, row.names = colnames(gene.mat))
			gene.mat <- trim.state.mat(gene.mat, gene.info, upstream.buffer, 
			downstream.buffer, start.feature, end.feature)
	
			#get distance matrix for chromatin states and SNPs
			d.chrom.mat <- get.hamming.mat(gene.mat)
			
			#get the SNPs for the same region
			gene.chr <- unique(gene.info[,3])
			
			if(gene.chr != "Y" && get.snps){
				snp.mat <- get.variants(chr = gene.chr, 
				start = min(as.numeric(rownames(gene.mat))), 
				end = max(as.numeric(rownames(gene.mat))), strains = var.strains)
			}else{
				snp.mat <- NULL	
			}
			
			if(length(snp.mat) > 0 && nrow(snp.mat) > 0){
				d.snp.mat <- get.hamming.mat(snp.mat[,var.strains])
				#test whether chromatin state is related to expression
				expr.state.comp <- compare.expr.chrom(d.chrom.mat, treat.expr)
				
				#test whether SNPs are related to expression
				expr.snp.comp <- compare.expr.chrom(d.snp.mat, treat.expr)
				
				if(nrow(snp.mat) > 0){
					snp.tree <- snps2tree(as.matrix(snp.mat[,var.strains]))
				}

			}else{
				snp.tree = NULL
				expr.state.comp <- NULL
				expr.snp.comp <- NULL
			}
					
			if(length(which(is.na(gene.mat))) == 0){
				chrom.tree <- snps2tree(gene.mat)
			}else{
					chrom.tree = NULL
			}
			
			result <- list(gene.info, gene.mat, treat.expr[med.order], snp.tree, chrom.tree, 
			expr.state.comp, expr.snp.comp)
			names(result) <- c("gene.info", "state.matrix", "control.expression", "snp.tree", 
			"chrom.state.tree", "expression.chromatin.cor.p", "expression.snp.cor.p")
		}else{
			result <- NULL
		}
	}else{
		result <- NULL	
	}
	
	
	return(result)
}
