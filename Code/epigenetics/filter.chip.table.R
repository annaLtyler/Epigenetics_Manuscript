filter.chip.table <- function(chip.table, FC = NULL, fdr = NULL){
	if(!is.null(FC)){chip.table <- chip.table[which(abs(chip.table[,"Fold"]) > FC),,drop=FALSE]}
	if(!is.null(fdr)){chip.table <- chip.table[which(abs(chip.table[,"FDR"]) <= fdr),,drop=FALSE]}		
	return(chip.table)
}