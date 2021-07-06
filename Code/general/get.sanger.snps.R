#this function retrieves sanger SNPs for a given location in the DO plus D2


get.sanger.snps <- function(chr, start.pos, end.pos){
	
	base.url = "www.sanger.ac.uk/sanger/Mouse_SnpViewer_Export/SNPs/rel-1505?sv=31;sn=8589410303;st=156274698;loc="
	location = paste0(chr, ":", start.pos, "-", end.pos)
	file.format = ";format=csv"
	full.url <- paste0(base.url, location, file.format)
	
	con <- curl(full.url)
	snp.table <- readLines(con)
	close(con)

	snp.colnames <- unlist(strsplit(snp.table[1], ","))

	snp.mat <- matrix(NA, ncol = length(snp.colnames), nrow = (length(snp.table)-1))
	colnames(snp.mat) <- snp.colnames
	
	for(i in 1:(length(snp.table)-1)){
		snp.mat[i,] <- unlist(strsplit(snp.table[(i+1)], ","))
		}
		
	return(snp.mat)
}