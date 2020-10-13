# This function calls doEdgeRAnalysis to compare expression
# values between groups (either treatment/control or strain
# comparisons) These comparisons are saved as qlf objects. 
# The function then calls EdgeR.to.table to generate human-
# readable tables from the analysis. These tables are also 
# filtered by pvalue and by fold change
ProcessAllRNA <- function(exp, pval = NULL, lFC = NULL){

	#set up directories
	all.dir <- dir.setup(exp)
	#assign the directory names to variables
	for(i in 1:length(all.dir)){assign(names(all.dir)[[i]], all.dir[[i]])} 

	#============================================================
	#read in RNA data and set up comparison groups
	#============================================================
	rna.base.dir <- strsplit(rna.data.dir, exp)[[1]]
	setwd(rna.base.dir)
	
	# get big file with all the strains
	strain.table=as.matrix(read.delim("5StrainsEffCts.txt", header=TRUE, sep="\t"))
	strain.col <- substr(colnames(strain.table), 1,2)
	strains <- sort(unique(strain.col))
	treat.col <- substr(colnames(strain.table), 3,3)
	treatments <- sort(unique(c(treat.col)))
	
	
	#============================================================
	#compare all strains within a treatment
	#============================================================
	
	if(exp == "strains"){
				
		strain.pairs <- pair.matrix(strains)
		for(s in 1:nrow(strain.pairs)){
			cat(strain.pairs[s,], "\n")
			strain1.locale <- which(strain.col == strain.pairs[s,1])
			strain2.locale <- which(strain.col == strain.pairs[s,2])
			for(tr in 1:length(treatments)){
				cat("\t", treatments[tr], "\n")
				treat.locale <- which(treat.col == treatments[tr])
				strain1.treat <- intersect(strain1.locale, treat.locale)
				strain2.treat <- intersect(strain2.locale, treat.locale)
				file.label <- paste(paste(strain.pairs[s,], collapse = ".v."), treatments[tr], sep = "_")
				count.table <- strain.table[,c(strain1.treat, strain2.treat)]
				groupV <- strain.col[c(strain1.treat, strain2.treat)]
				qlf.obj <- doEdgeRAnalysis(count.table, groups = groupV)
				
				setwd(rna.data.dir)
				saveRDS(qlf.obj, paste(file.label, ".RData", sep = ""))
				
				setwd(rna.results.dir)
				#Collect additional information on the diff exp genes
				#and write to table
				EdgeR.to.table(qlf.obj, result.name = paste("Diff.Exp.Genes.", file.label, ".txt", sep = ""), pval = pval, lFC = lFC)
				}
			}
		} #end processing of strains
	
	
	
	#===============================================================
	# Specify directories for the treatments
	#===============================================================
	if(exp == "TvC"){
	
	
	for(s in 1:length(strains)){
		cat(strains[s], "\n")
		strain.locale <- which(strain.col == strains[s])
		treat1.locale <- which(treat.col == treatments[1])
		treat2.locale <- which(treat.col == treatments[2])
		
		strain.treat1 <- intersect(strain.locale, treat1.locale)
		strain.treat2 <- intersect(strain.locale, treat2.locale)
		
		file.label <- paste(strains[s], paste(treatments, collapse = ".v."), sep = "_")
		count.table <- strain.table[,c(strain.treat1, strain.treat2)]
		groupV <- treat.col[c(strain.treat1, strain.treat2)]
		
		#test for differentially expressed genes using EdgeR
		qlf.obj <- doEdgeRAnalysis(count.table, groups = groupV)
		setwd(rna.data.dir)
		saveRDS(qlf.obj, paste(file.label, ".RData", sep = ""))
	
		#Collect additional information on the diff exp genes
		#and write to table
		setwd(rna.results.dir)
		EdgeR.to.table(qlf.obj, result.name = paste("Diff.Exp.Genes.", file.label, ".txt", sep = ""), pval = pval, lFC = lFC)
		
		} #end looping through strains
	} #end if exp == TvC	
} #end function
