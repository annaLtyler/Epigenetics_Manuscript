# chr.data <- read.table("~/Documents/Epigenetics/Data/genes/gene_coords.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)


get.gene.chr <- function(gene.list, chr.data, max.chr = 19){

	require(usefulScripts)
	
	#find the chromosome of each gene in the query
	query.table <- matrix(NA, ncol = 1, nrow = length(gene.list))
	rownames(query.table) <- gene.list
	for(g in 1:length(gene.list)){
		gene.locale <- which(chr.data[,3] == gene.list[g])
		if(length(gene.locale) > 0){
			query.table[g,1] <- chr.data[gene.locale,5]
			}
		}

	query.table[which(query.table[,1] == "X"),] <- max.chr + 1
	query.table[which(query.table[,1] == "Y"),] <- max.chr + 2
	query.table[which(query.table[,1] == "MT"),] <- max.chr + 3
	
	chr.order <- order(as.numeric(query.table[,1]))

	query.table[which(query.table[,1] == max.chr + 1),] <- "X"
	query.table[which(query.table[,1] == max.chr + 2),] <- "Y"
	query.table[which(query.table[,1] == max.chr + 3),] <- "MT"

	sorted.table <- as.matrix(query.table[chr.order,])
	return(sorted.table)
	}