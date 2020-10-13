#This function takes in an enrichment table returned from gProfileR
#and *write* the results to a file


write.enrichment <- function(enrichment, num.rows = 10, text.size = 1.5, filename = "Enrichment.txt"){

	if(nrow(enrichment) == 0){
		message("No Significant Enrichment")
		return()
		}
		
	num.rows <- min(c(num.rows, nrow(enrichment)))
	columns.to.write <- c("term.name", "term.size", "query.size", "overlap.size", "p.value", "intersection")

	sub.table <- enrichment[1:num.rows,columns.to.write]
	colnames(sub.table) <- c("term", "term.size", "query.size", "overlap", "p.value", "genes")

	write.table(sub.table, filename, quote = FALSE, sep = "\t", row.names = FALSE)

	}	

