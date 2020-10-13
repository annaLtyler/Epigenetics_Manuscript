ProcessChIP <- function (mark.type, group1, group2, analysis, plot.results = FALSE){
	
	contrast.table <- dba.show(analysis, bContrasts = TRUE)
	direction.flag <- 0
	
	for(i in 1:length(group1)){
		cat("\t", group1[i], "vs.", group2[i], "\n")
		#find the location of the comparison in the contrast table
		contrast.idx <- intersect(which(contrast.table[,"Group1"] == group1[i]), which(contrast.table[,"Group2"] == group2[i]))
		
		if(length(contrast.idx) == 0){
			warning("contrast in wrong direction")
			contrast.idx <- intersect(which(contrast.table[,"Group2"] == group1[i]), which(contrast.table[,"Group1"] == group2[i]))
			direction.flag <- 1
			}

		## Return GRanges with FDR < 0.1
		chip.GR <- dba.report(analysis, contrast=contrast.idx)
	  
		#make results plots	
		if(direction.flag == 0){
			plot.label <- paste0(mark.type, "_", group1[i], "_", group2[i])
			}else{
			plot.label <- paste0(mark.type, "_", group2[i], "_", group1[i])	
			}
			plot.label <- gsub(":", "_", plot.label)
			
			
			#Make dataframe from GRanges
			chip.df <- data.frame(chr = seqnames(chip.GR),start = start(chip.GR)-1,end = end(chip.GR), mcols(chip.GR))
			## Take out Y chromosome and MT chromosome
			chip.df = chip.df[chip.df$chr != "Y",]
			chip.df = chip.df[chip.df$chr != "MT",]


			if(plot.results){
			pdf(file = paste0(plot.label, ".pdf"))
			hist(chip.GR$p.value, main = paste("Histogram of", plot.label, "P-values"), breaks = 100)
		  	
		
			dba.plotMA(analysis, contrast = contrast.idx, bXY = TRUE)
			legend("topleft", col = "maroon1", legend = "Significant", pch = 16, cex = 0.7)
			
			dba.plotPCA(analysis, contrast = contrast.idx)
			
			pvals <- dba.plotBox(analysis, contrast = contrast.idx)
	
			# dba.plotHeatmap(analysis)
	
			dba.plotHeatmap(analysis, contrast = contrast.idx, correlations = FALSE, scale = "row")
			dev.off()
			}
		
		
		file.name <- paste0(plot.label, ".txt")
		
		#if the analysis is in the opposite direction as the RNA
		#switch the direction
		if(direction.flag == 1){
			colnames(chip.df)[5:6] <- colnames(chip.df)[6:5]
			chip.df[,5:6] <- chip.df[,6:5]
			}else{
			#The ChIP fold change is negative if the concentration increases from group1 to group2	
			chip.df[,7] <- chip.df[,7]*-1
			}

			
		write.table(chip.df, file = file.name, row.names = FALSE, sep = "\t", quote = FALSE)
	
		} #end looping through group comparisons

} #end function

