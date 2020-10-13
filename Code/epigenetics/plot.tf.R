#This function takes in an enrichment table returned from gProfileR
#and *plots* the results in an easier-to-read table


plot.tf <- function(tf.table, num.rows = 10, text.size = 1.5, plot.label = "TF Enrichment"){

	if(length(tf.table) == 0){
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, 1))
		text(x = 0.5, y = 0.5, "No Significant TFs")
		return()
		}
		
	num.rows <- min(c(num.rows, nrow(tf.table)))

	sub.table <- tf.table[1:num.rows,]
	sub.table[,"p.value"] <- signif(sub.table[,"p.value"], 2)
	colnames(sub.table) <- gsub("Genes", "G", colnames(sub.table))
	colnames(sub.table) <- gsub("Genome", "Gen", colnames(sub.table))

	x.pts <- c(0.3, segment.region(0.45, 1, (ncol(sub.table) - 1), alignment = "ends"))
	if(nrow(sub.table) > 1){
		y.pts <- segment.region(1, 0, nrow(sub.table), alignment = "ends")
		}else{
		y.pts <- 0.5	
		}

	plot.new()
	plot.window(xlim = c(0, 1), ylim = c(0, 1))
	#write the column and row names of the matrix
	par(xpd = TRUE)
	text(x = 0.5, y = 1.1, plot.label)
	text(x = x.pts, y = rep(y.pts[1], (length(x.pts)-1)), colnames(sub.table), adj = 1, cex = text.size)
	for(d in 1:dim(enrichment)[1]){
		text(x = x.pts, y = rep(y.pts[(d+1)], (length(x.pts)-1)), sub.table[d,], col = "black", adj = 1, cex = text.size)
		}
	par(xpd = FALSE)
	}	

