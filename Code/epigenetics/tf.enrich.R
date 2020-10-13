#Get TRANSFAC enrichment for predicted TF binding sites
#use GATHER website: http://gather.genome.duke.edu
	# 1	Use your favorite programming language (or wget) and send GET or POST requests to:
	# 2	http://gather.genome.duke.edu/?cmd=report&gene_box=<genes>&tax_id=<organism>&annot_type=<annot_type>&network=<network>&homologs=<homologs>
	# 3	Parameters:
# cmd: Must be "report".
# gene_box: Gene IDs separated by spaces or commas.
# tax_id: Entrez Taxonomy ID for your organism, e.g. 9606 for human.
# annot_type: Specifies the type of annotation returned. Must be one of: gene_ontology for Gene Ontology annotations, words for MEDLINE Words, mesh for MeSH, kegg for KEGG Pathway, proteins for Protein Binding, literature for Literature Net, mirna for miRNA, transfac for TRANSFAC, maploc for Chromosome.
# network: Should be 0 or 1 indicating whether annotations from a literature network should be used.
# homologs: Should be 0 or 1 indicating whether annotations from homologs should be used.
	# 4	e.g. http://gather.genome.duke.edu/?cmd=report&gene_box=e2f1+e2f3+myc&tax_id=9606&annot_type=gene_ontology&network=0&homologs=1 The server will return the results as a tab-delimited text file. Please add a 2s delay between queries so that the server is not saturated with requests. We do monitor for excessive use and will throttle or cut off a connection if it is degrading the quality of service for others. If you need to run a large batch of queries, contact us.
#gene.list <- c("Treh", "Cyp4a12a", "Serpina1e", "Hsd3b5", "rpina4−ps1", "Cyp2d9", "Mup21", "Slc22a28", "Derl3", "Slco1a1", "Sult1e1", "C6", "Serpine2", "Arsa", "Acsm2", "Usp2", "C8a", "C8b", "Cd9", "Cyp4a31", "Mug1", "Slc17a8", "Cyp4a14", "Srd5a1", "Hsd3b3", "Pdilt", "Alas2", "Crlf2", "Cyp3a11", "Ugt2b5", "C9", "Ociad2", "Il1rap", "Slc13a2", "Extl1", "Tspan33", "Ugt2b38", "Acad10", "Mycl", "Cyp3a41a", "Gas6", "Col15a1", "Slco1a4", "Frmd4b", "Ddx60", "Ifit3", "Abcc4", "Cyp17a1", "Cyp2b9", "Cyp2c37", "Arntl", "Cmpk2", "Gbp2", "Trim24", "Ifit1", "Serpinb1a", "Hamp2", "Arrdc4", "Mme", "Gbp3", "Gbp6", "Rcan2", "Fam46c", "Rbp1", "Lpin1", "Cd74", "Tpm2", "Ly6a", "H2−Eb1", "Cyp3a44", "Slc22a27", "Hao2", "Nipal1", "Atp6v0d2", "Cyp2b13", "Sult2a2", "Sult2a1", "Col4a3", "A1bg", "Cyp2a22", "Slc22a26", "Sult2a5", "Cyp3a41b", "Cyp3a1")

tf.enrich <- function(gene.list, dest.dir = ".", filename, sig.level = 0.05){

	if(dest.dir == "."){dest.dir <- getwd()}
	file.dest <- paste(dest.dir, filename, sep = "/")

	gene_box = paste(gene.list, collapse = ",")	
	query.line <- paste("http://gather.genome.duke.edu/?cmd=report&gene_box=", gene_box, "&tax_id=10090&annot_type=transfac&network=0&homologs=1", sep = "")
	
	download.file(query.line, file.dest)

	
	transfac.data <- try(read.table(filename, sep = "\t", header = FALSE, stringsAsFactors = FALSE), silent = TRUE)
	
	if(class(transfac.data) == "try-error"){return()}

		p.col = 9
		name.col = 2
	
		pval <- exp(transfac.data[,p.col]*-1)
		sig.locale <- which(pval <= sig.level)
		
		sig.tf <- cbind(transfac.data[sig.locale,2:7], pval[sig.locale])
		colnames(sig.tf) <- c("Annotation", "Genes.Ann", "QGenes.Ann", "QGenes.NoAnn", "Genome.Ann", "Genome.NoAnn", "p.value")	
		
		return(sig.tf)	
}









