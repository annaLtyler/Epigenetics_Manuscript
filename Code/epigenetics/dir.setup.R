dir.setup <- function(exp = c("TvC", "strains")){

	#if exp isn't changed, default to strains
	if(length(exp) == 2){
		exp <- "strains"
		}

	base.dir <- "~/Documents/Epigenetics"
	
	rna.data.dir <- paste(base.dir, "Data/RNASeq", exp, sep = "/")
	chip.data.dir <- paste(base.dir, "Data/ChIP", exp, sep = "/")
	
	rna.results.dir <- paste(base.dir, "Results/RNASeq", exp, sep = "/")
	chip.results.dir <- paste(base.dir, "Results/ChIP", exp, sep = "/")
	
	combined.results.dir <- paste(base.dir, "Results/combined", exp, sep = "/")

	dir.list <- list("base.dir" = base.dir, "rna.data.dir" = rna.data.dir, "chip.data.dir" = chip.data.dir, "rna.results.dir" = rna.results.dir, "chip.results.dir" = chip.results.dir, "combined.results.dir" = combined.results.dir)

	return(dir.list)

	}