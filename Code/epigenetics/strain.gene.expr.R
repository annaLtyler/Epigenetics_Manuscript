# rna.table <- as.matrix(read.table("~/Documents/Epigenetics/Data/RNASeq/9StrainsEffCts.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t"))
# rna.table <- log(round(rna.table))
# rna.table[which(!is.finite(rna.table))] <- 0
# mus = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

strain.gene.expr <- function(rna.table, gene.name, comp = c("C", "C.v.T"), mart){

	require(biomaRt)
	require(usefulScripts)
	
	gene.id <- getBM("ensembl_gene_id", "external_gene_name", gene.name, mus)
	
	
	strains <- substr(colnames(rna.table), 1, 2)
	u_strain <- sort(unique(strains))
	treat <- substr(colnames(rna.table), 3, 3)
	u_treat <- sort(unique(treat))
	
	gene.idx <- which(rownames(rna.table) == gene.id[1,1])
	
	if(length(gene.idx) == 0){
		stop("Couldn't find gene.")
		}
	
	get.strain.data <- function(strain, treatment){
		strain.locale <- which(strains == strain)
		data.locale <- intersect(strain.locale, which(treat == treatment))
		return(rna.table[gene.idx, data.locale])
		}

	ctrl.list <- lapply(u_strain, function(x) get.strain.data(x, "C"))
	names(ctrl.list) <- u_strain
	
	if(length(comp) > 1 || comp == "C"){
		plot.order <- order(unlist(lapply(ctrl.list, median)), decreasing = FALSE)
		boxplot(ctrl.list[plot.order], ylab = "log(TPM)", xlab = "Strain")
		}else{
		treat.list <- lapply(u_strain, function(x) get.strain.data(x, "T"))
		names(treat.list) <- u_strain
		
		no.vals <- which(lapply(treat.list, length) == 0)
		ctrl.list <- ctrl.list[-no.vals]
		treat.list <- treat.list[-no.vals]
		
		plot.order <- order(unlist(lapply(ctrl.list, median)), decreasing = FALSE)
		plot.grouped.boxes(group.list = list(ctrl.list[plot.order], treat.list[plot.order]), group.labels = c("Control", "Treatment"), main = gene.name)
		}

	
	
}