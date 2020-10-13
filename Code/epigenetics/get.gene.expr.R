#This function gets the expression of a gene
#and finds averages for different classes
#(strains or treatment) if desired


get.gene.expr <- function(rnaseq.gene.info, rna.seq, gene.name, strain.name, treatment.name, 
col.table, average.replicates = TRUE){
	
	strains <- substr(colnames(rna.seq), 1, 2)
	treatments <- substr(colnames(rna.seq), 3,3)
	
	name.locale <- which(rnaseq.gene.info[,"external_gene_name"] == gene.name)[1]
	gene.id <- rnaseq.gene.info[name.locale,"ensembl_gene_id"]
	gene.locale <- which(rownames(rna.seq) == gene.id)
	
	if(length(strain.name) == 1){
		strain.locale <- which(strains == strain.name)
	}else{
		strain.locale <- which(!is.na(order.strains(strains, strain.name, col.table)))
	}
	treatment.locale <- which(treatments %in% treatment.name)
	id.locale <- intersect(strain.locale, treatment.locale)
	
	if(length(id.locale) > 0){
		expr <- rna.seq[gene.locale,id.locale]
		non.rep.id <- substr(names(expr), 1,3)
		u_ids <- unique(non.rep.id)
		if(average.replicates){
			result <- rep(NA, length(u_ids))
			names(result) <- u_ids
			}else{
				result <- matrix(NA, ncol = length(u_ids), nrow = length(non.rep.id)/length(u_ids))
				colnames(result) <- u_ids	
				if(length(strain.name) == 1){ #if we are only looking at one strain, keep track of the ID
					rownames(result) <- names(expr)
				}
			}
		for(i in 1:length(u_ids)){
			id.locale <- which(non.rep.id == u_ids[i])
			if(average.replicates){
				result[i] <- mean(expr[id.locale])
				}else{
				result[,i] <- expr[id.locale]
				}
			}
		}
	
	return(result)
	
}