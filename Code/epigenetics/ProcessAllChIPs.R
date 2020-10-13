#This function calls ProcessChIP to generate data 
#tables with gene information for genes that have 
#differential histone marks. These tables can then
#be read in by other functions for further processing.
#exp can either be TvC or strains

ProcessAllChIPs <- function(exp = "TvC"){	
	
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
	treatments <- c("control", "treated")
	
	if(exp == "TvC"){
		#treatment vs. control
		group1  <- paste(strains, treatments[1], sep = ":")
		group2 <- paste(strains, treatments[2], sep = ":")
		}else{			
		#individual strains
		strain.pairs <- pair.matrix(strains)
		group1  <- strain.pairs[,1]
		group2 <- strain.pairs[,2]
		}
			
		#get the names of the files. Each one looks at a 
		#different epigenetic marker in treatment v. control
		#for the strains listed above. We want to compare
		#treatment v. control for each strain and each mark.
		 
		setwd(chip.data.dir)
		chip.names = list.files()
			
		for(ch in 1:length(chip.names)){
			setwd(chip.data.dir)
			mark.type <- strsplit(chip.names[ch], "_")[[1]][2]
			cat(mark.type, "\n")
		
		 	#This loads the variable "analysis." 
		 	#It is a DBA object from a DiffBind analysis.
		  	load(chip.names[ch])
		
			results.dir <- paste(chip.results.dir, mark.type, sep = "/")
			if(!file.exists(results.dir)){system(paste("mkdir", results.dir))}
			setwd(results.dir)
			
			## ChIP processing
			ProcessChIP(mark.type, group1, group2, analysis)
		  	}
	 }
	