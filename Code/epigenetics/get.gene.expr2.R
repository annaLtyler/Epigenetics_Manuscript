

get.gene.expr2 <- function(gene.name, gene.info.table, rna.seq, col.table, C.or.T = "C"){
	
	gene.locale <- which(gene.info.table[,2] == gene.name)
	gene.id <- gene.info.table[gene.locale,1]
	id.locale <- which(rownames(rna.seq) == gene.id)
	
	exp.locale <- grep("T", colnames(rna.seq))
	
	if(C.or.T == "C"){
		exp.locale <- setdiff(1:ncol(rna.seq), exp.locale)
		}
	
	exp.rna.seq <- rna.seq[,exp.locale]
	strains <- substr(colnames(exp.rna.seq), 1,2)
	u_strains <- unique(strains)
	
	strain.expr <- lapply(u_strains, function(x) log(exp.rna.seq[id.locale,which(strains == x)]))
	names(strain.expr) <- u_strains
	return(strain.expr)
	
}