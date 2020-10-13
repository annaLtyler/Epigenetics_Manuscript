#This function runs the overlap analysis for the chip and rna data
#by calling CorrelateRNAwChIP
#when looking for histone marks near differentially expressed genes, the upstream.bp and downstream.bp indicate how far outside of the gene start and end to look for epigenetic marks associated with that gene.
# exp = "TvC" or "strains"
# rna.lFC = minimum log fold change for RNA
# rna.pval = maximum p value for RNA
# chip.FC = minimum fold change for chip
# chip.fdr = maximum fdr for chip (fdr limit of 0.1 was used in DiffBind Analysis, so an fdr of this value will not filter the chip)


Combine_RNA_ChIP <- function(exp = "TvC", upstream.bp = 5000, downstream.bp = 5000, rna.lFC = NULL, rna.pval = NULL, chip.FC = NULL, chip.fdr = NULL){

		
	#===============================================================
	# Specify directories
	#===============================================================
	#get all the directories from dir.setup.R
	all.dir <- dir.setup(exp)
	#assign the directory names to variables
	for(i in 1:length(all.dir)){assign(names(all.dir)[[i]], all.dir[[i]])} 
	
		
	#===============================================================
	#define the groups to be compared
	#===============================================================
	strains <- c("B6", "CAST", "DBA", "PWK", "WSB")
	strain.nicknames <- c("B6", "CA", "DB", "PK", "WS")
	treatments <- c("control", "treated")
	treatment.nicknames <- c("C", "T")
	mark.types <- c("k27", "me1", "me3")


	if(exp == "TvC"){
		#===============================================================
		# Do analysis of ChIP data, treatment vs. control
		#===============================================================
		
		#===============================================================
		# Compare the RNASeq and ChIP data
		# for each strain and each epigenetic mark
		#===============================================================
	
		 for(s in 1:length(strains)){
			cat(strains[s], "\n")
			#get the RNA data
			setwd(rna.results.dir)
			rna.file <- get.files(want = strain.nicknames[s], dont.want = "pdf")
			rna.data <- read.table(rna.file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
		 	rna.data <- filter.rna.table(rna.data, lFC = rna.lFC, pval = rna.pval)
		 	
		 	for(m in 1:length(mark.types)){
		 		cat("\t", mark.types[m], "\n")
				results.file <- paste(mark.types[m], "_", strains[s], ".txt", sep = "")
				
		 		#get the CHiP data
			 	setwd(paste(chip.results.dir, mark.types[m], sep = "/"))
		 		chip.file <- get.files(want = c(strains[s], ".txt"))
		 		chip.data <- read.table(chip.file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
		 		chip.data <- filter.chip.table(chip.data, FC = chip.FC, fdr = chip.fdr)
		 		
				setwd(paste(combined.results.dir, mark.types[m], sep = "/"))
				
				##Correlate RNA with ChIP
				ChIP_RNA_Overlap(rna.data, chip.data, outfile = results.file, upstream.bp, downstream.bp)
		 		} #end looping through epigenetic marks
		 	} #end looping through strains
		 
		  
		#===============================================================
		# Finish Treatment v. Control
		#===============================================================	
		}
	
	#===============================================================
	# Compare strains (control data only)
	#===============================================================

	if(exp == "strains"){
				
		#===============================================================
		# Process the ChIP files
		#===============================================================
			
		#define the groups as they were defined in the RNASeq data
		setwd(rna.results.dir)
		filenames <- get.files(dont.want = ".pdf")
	
		group1.nicknames <- substr(filenames, 16, 17)
		group2.nicknames <- substr(filenames, 21, 22)
		
		#change the names in the groups to the full names
		#these are used in the ChIP analyses
		group1 <- group1.nicknames
		group2 <- group2.nicknames
		for(i in 1:length(strains)){
			group1[which(group1 == strain.nicknames[i])] <- strains[i]
			group2[which(group2 == strain.nicknames[i])] <- strains[i]
			}
		
		#===============================================================
		# Compare the RNASeq and CHiP data
		# for each strain comparison and each 
		# epigenetic mark
		#===============================================================
			
		rna.comparisons <- apply(cbind(group1.nicknames, group2.nicknames), 1, function(x) paste(x, collapse = ".v."))
		chip.comparisons <- apply(cbind(group1, group2), 1, function(x) paste(x, collapse = "_"))
			
		 for(s in 1:length(rna.comparisons)){
			cat(rna.comparisons[s], "\n")

			#get the RNA data
			setwd(rna.results.dir)
			rna.file <- get.files(want = c(rna.comparisons[s], "_C"), dont.want = ".pdf")
			rna.data <- read.table(rna.file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	 	
		 	for(m in 1:length(mark.types)){
		 		cat("\t", mark.types[m], "\n")
				results.file <- paste(mark.types[m], "_", chip.comparisons[s], ".txt", sep = "")
				
		 		#get the CHiP data
			 	setwd(paste(chip.results.dir, mark.types[m], sep = "/"))
		 		chip.file <- get.files(want = c(chip.comparisons[s], ".txt"))
		 		chip.data <- read.table(chip.file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
		 		
				setwd(paste(combined.results.dir, mark.types[m], sep = "/"))
			
				##Correlate RNA with ChIP
				ChIP_RNA_Overlap(rna.data, chip.data, outfile = results.file, upstream.bp, downstream.bp)
		 		} #end looping through epigenetic marks
		 	} #end looping through strains		 
		 } #end case for exp == "strains"
	} #end function
	
