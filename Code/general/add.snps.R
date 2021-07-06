#This function adds SNPs to a sequence based on a Sanger
#SNP csv output and a sequence from b6.sequence()


add.snps <- function(b6.seq, sanger.snps){
	require("stringr")
	
	csq.col <- grep("Csq", colnames(sanger.snps))
	no.csq <- sanger.snps[,-csq.col]
	
	just.snps <- no.csq[,6:ncol(no.csq),drop=FALSE]
	strains <- gsub("X", "", colnames(just.snps))
	
	#make a snp table for each strain
	seq.table <- matrix(b6.seq[,2], ncol = length(strains), nrow = nrow(b6.seq), byrow = FALSE)
	for(i in 1:length(strains)){
		snp.locale <- which(just.snps[,i] != "-")
		snp.coord <- sanger.snps[snp.locale,2]
		strain.snps <- str_to_lower(as.vector(just.snps[snp.locale,i]))
		change.locale <- match(snp.coord, b6.seq[,1])
		seq.table[change.locale,i] <- strain.snps
		}
	
	colnames(seq.table) <- strains
	colnames(b6.seq) <- c("position", "B6")
	final.table <- cbind(b6.seq, seq.table)

	#remove symbols that aren't allowed
	just.nuc <- final.table[,2:ncol(final.table)]
	not.allowed <- which(unlist(apply(just.nuc, 1, function(x) length(which(!x %in% c("a", "c", "t", "g"))))) > 0)
	if(length(not.allowed) > 0){
		final.table <- final.table[-not.allowed,]
		}

	return(final.table)
	
}