#This function returns the position of a gene in rnaseq.gene.info
#from epigen.states.and.expr

get.gene.i <- function(rnaseq.gene.info, gene.name = "Akr1c18"){
	name.locale <- which(rnaseq.gene.info[,"external_gene_name"] == gene.name)
	gene.id <- unique(rnaseq.gene.info[name.locale,"ensembl_gene_id"])
	id.locale <- which(rownames(rna.seq) == gene.id)
	return(id.locale)
	}
	
