filter.rna.table <- function(rna.table, lFC = NULL, pval = NULL){
	if(!is.null(lFC)){rna.table <- rna.table[which(abs(rna.table[,"logFC"]) > lFC),,drop=FALSE]}
	if(!is.null(pval)){rna.table <- rna.table[which(abs(rna.table[,"PValue"]) <= pval),,drop=FALSE]}		
	return(rna.table)
}