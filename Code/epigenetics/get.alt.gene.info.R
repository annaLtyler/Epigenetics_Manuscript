#This function retrieves gene information from rnaseq.gene.info
#from epigen.states.and.expr
#it uses a region other than the gene start and stop, though,
#so you can look a putative regulatory regions compared with
#gene expression

get.alt.gene.info <- function(rnaseq.gene.info, rna.seq, gene.id, start.position, end.position, upstream.buffer, downstream.buffer){

	var.strains <- c("129S1/SvImJ", "A/J", "C57BL/6J", " CAST/EiJ", "DBA/2J", "NOD/ShiLtJ", "NZO/HlLtJ", "PWK/PhJ", "WSB/EiJ")
	
	gene.locale <- which(rnaseq.gene.info[,1] == gene.id)
	gene.expr.locale <- which(rownames(rna.seq) == gene.id)
	
	if(length(gene.locale) > 0){
		gene.info <- rnaseq.gene.info[gene.locale,,drop=FALSE]
	
		#use larger buffers then trim to make sure we get states
		strain.states <- lapply(all.bed, function(x) get.states.spec.pos(gene.info, x, upstream.buffer = 50000, downstream.buffer = 50000, start.position, end.position))
		names(strain.states) <- strain.names
			
		if(sum(unlist(lapply(strain.states, nrow))) > 0){

			#pull out expression values for controls
			gene.expr <- rna.seq[gene.expr.locale,]	
			labels <- substr(names(gene.expr), 1, 3)
			u_labels <- sort(unique(labels))
			binned.expr <- lapply(u_labels, function(x) as.vector(gene.expr[which(labels == x)]))
			names(binned.expr) <- u_labels
			c.locale <- which(substr(names(binned.expr), 3, 3) == "C")
			c.expr <- binned.expr[c.locale]
			med.order <- order(unlist(lapply(binned.expr[c.locale], median)))

			
			#line up the states for the different strains and trim to requested region
			gene.mat <- align.states(strain.states)
			gene.mat <- trim.state.mat(gene.mat, gene.info, upstream.buffer, downstream.buffer, start.feature, end.feature)
	
			#get distance matrix for chromatin states and SNPs
			d.chrom.mat <- get.hamming.mat(gene.mat)
			
			#get the SNPs for the same region
			gene.chr <- unique(gene.info[,3])
			if(gene.chr != "Y"){
			snp.mat <- get.variants(chr = gene.chr, start = min(as.numeric(rownames(gene.mat))), end = max(as.numeric(rownames(gene.mat))), strains = var.strains)
			}else{
			snp.mat <- NULL	
			}
			
			if(length(snp.mat) > 0 && nrow(snp.mat) > 0){
				d.snp.mat <- get.hamming.mat(snp.mat[,var.strains])
				#test whether chromatin state is related to expression
				expr.state.comp <- compare.expr.chrom(d.chrom.mat, c.expr)
				
				#test whether SNPs are related to expression
				expr.snp.comp <- compare.expr.chrom(d.snp.mat, c.expr)
				
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
			
			result <- list(gene.info, gene.mat, c.expr[med.order], snp.tree, chrom.tree, expr.state.comp, expr.snp.comp)
			names(result) <- c("gene.info", "state.matrix", "control.expression", "snp.tree", "chrom.state.tree", "expression.chromatin.cor.p", "expression.snp.cor.p")
			}
		}else{
		result <- NULL	
		}
	
	
	return(result)
}
