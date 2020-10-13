#This function finds the overlap between genes that are differentially
#expressed and genes that have epigenetic modifications. The overlap with
#epigenetic modifications can be restricted to just the coding region, or 
#up and downstream regions from the coding region using upstream.bp and
#downstream.bp.

ChIP_RNA_Overlap <- function(rna.data, chip.data, outfile, upstream.bp = 5000, downstream.bp = 5000){

	## Convert text file to GRanges
	chip.GR = GRanges(seqnames = chip.data$chr,
			ranges = IRanges(start = chip.data$start,
			end = chip.data$end),
			mcols = chip.data[,5:9])


	# Convert to GRanges 
	RNA.GR = GRanges(seqnames = rna.data$Chr,
			ranges = IRanges(start = rna.data$Start - upstream.bp,
			end = rna.data$End + downstream.bp), mcols = rna.data[,c(1:2, 6:9)])


	##Intersect intervals with GRanges.
	ol = findOverlaps(query = RNA.GR, subject = chip.GR)

	if(length(ol) == 0){return()}

	#subset RNA GRanges obj to get regions that overlap with chip data
	ol.rna = RNA.GR[queryHits(ol),]
	
	#check to see how the overlaps match with the chip ranges
	ol.chip <- chip.GR[subjectHits(ol),]
	
	ol.rna <- as.data.frame(ol.rna)
	colnames(ol.rna) <- c("Chr", "start", "end", "width", "strand", "EnsemblID", "Gene", "RNA_logFC", "RNA_logCPM", "RNA_F", "RNA_PValue")
	
	ol.chip <- as.data.frame(ol.chip)
	colnames(ol.chip) <- c("Chr", "start", "end", "width", "strand", "Conc_B6.control", "Conc_B6.treated", "ChIP_logFC", "ChIP_p.value", "ChIP_FDR")
	
	full.table <- cbind(ol.rna[,c(1:3,6:9)], ol.chip[,c(6:10)])
	
	## Order by RNA logFC
	full.table = full.table[order(full.table$RNA_logFC),]

	## Add results metrics of # genes, # peaks, # overlaps
	filler <- rep("-", (nrow(full.table)-1))
	full.table$NumGenes = c(length(RNA.GR), filler)
	full.table$NumPeaks = c(length(chip.GR), filler)
	full.table$NumOverlaps = c(length(ol), filler)

	if(length(ol) > 5){
		correlation <- cor.test(ol.rna$RNA_logFC, ol.chip$ChIP_logFC, method = "pearson")
		model <- lm(ol.chip$ChIP_logFC~ol.rna$RNA_logFC)
		full.table$PearsonCor <- c(signif(correlation$estimate, 3), filler)
		full.table$PearsonCorP <- c(signif(correlation$p.value, 3), filler)
		pdf(paste0(strsplit(outfile, ".txt")[[1]], ".pdf"))
		plot(full.table$RNA_logFC, full.table$ChIP_logFC, xlab = "log FC RNA", ylab = "log FC ChIP")
		abline(h = 0, lty = 2, col = "gray")
		abline(v = 0, lty = 2, col = "gray")
		abline(model)
		dev.off()
		}else{
		full.table$PearsonCor <- c(NA, filler)
		full.table$PearsonCorP <- c(NA, filler)
		}
		

	write.table(full.table, row.names = FALSE, file = outfile, sep = "\t", quote = FALSE)
	
} 
	
	
