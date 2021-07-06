#This function downloads B6 sequence from ensembl

b6.sequence <- function(chr, start.pos, end.pos){
	
	base.url <- "http://genome.ucsc.edu/cgi-bin/das/mm10/dna?segment="
	region.tag <- paste0("chr", chr, ":", start.pos, ",", end.pos)
	full.url <- paste0(base.url, region.tag)
	
	download.file(full.url, "region.txt")
	b6.seq <- as.matrix(read.table("region.txt", fill = TRUE, stringsAsFactors = FALSE))
	seq.data <- b6.seq[,1]	
	no.data.rows <- grep("<", seq.data)
	data.rows <- setdiff(1:length(seq.data), no.data.rows)
	nucleotide <- unlist(strsplit(paste0(as.vector(seq.data[data.rows]), collapse = ""), ""))
	position <- start.pos:end.pos

	coord.table <- cbind(position, nucleotide)
	
	unlink("region.txt")
	
	return(coord.table)
}