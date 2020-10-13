get.methyl.by.gene <- function(gene.name, gene.info.table, methyl.data, 
methyl.id.table, upstream.buffer = 0, downstream.buffer = 0, strains = NULL, 
treatment = c("C", "T"), strain.means = TRUE){
	
	gene.locale <- which(gene.info.table[,2] == gene.name)
	if(length(gene.locale) == 0){
		return("No data")
		}
	
	treatment <- treatment[1]

	# gene.info.table[gene.locale,]
	gene.start <- gene.info.table[gene.locale[1],"start_position"] - upstream.buffer
	gene.end <- gene.info.table[gene.locale[1],"end_position"] + downstream.buffer
	gene.chr <- gene.info.table[gene.locale[1],"chromosome_name"]
	gene.strand <- as.numeric(gene.info.table[gene.locale[1],"strand"])

	if(is.null(strains)){strains <- unique(methyl.id.table[,"strain"])}

	methyl.table <- get.methyl.pos(all.methyl.data = methyl.data, 
	strain.name = strains, treatment.name = treatment, id.table = methyl.id.table, 
	chr = gene.chr, start.pos = gene.start - upstream.buffer, 
	stop.pos = gene.end + downstream.buffer, average.replicates = strain.means)

	return(methyl.table)

	}