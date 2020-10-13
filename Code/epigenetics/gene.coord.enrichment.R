# gene.coord <- read.table("~/Documents/Epigenetics/Data/genes/gene_coords.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)


gene.coord.enrichment <- function(gene.list, gene.coord, n.samples = 100, sig.val = 0.05, plot.label = "Significantly Enriched Locations"){

	library(usefulScripts)
		
	#find the number of genes per chromosome
	all.chr <- c(1:19, "X")
	gene.table <- matrix(NA, ncol = 1, nrow = length(all.chr))
	rownames(gene.table) <- all.chr
	colnames(gene.table) <- "num.genes"
	
	
	for(ch in 1:length(all.chr)){
		chr.locale <- which(gene.coord[,5] == all.chr[ch])
		gene.table[ch,1] <- length(unique(gene.coord[chr.locale,3]))
		}
	
	probs <- gene.table[,1]/sum(gene.table[,1])
	
	#find the chromosome of each gene in the query
	query.table <- get.gene.chr(gene.list, gene.coord, max.chr = 19)
	
	gene.dist <- table(query.table[,1])
	gene.counts <- gene.table; gene.counts[,1] <- 0
	gene.counts[names(gene.dist),] <- gene.dist
	
	expected.vals <- sum(gene.counts)*probs
	full.table <- cbind(gene.counts, expected.vals)
	
	
	#look for enrichment of location in the query list
	null.dist <- matrix(0, nrow = length(all.chr), ncol = n.samples)
	rownames(null.dist) <- all.chr
	
	for(n in 1:n.samples){
		loc.sample <- sample(all.chr, size = length(gene.list), prob = probs, replace = TRUE)
		counts <- table(loc.sample)
		null.dist[names(counts),n] <- counts
		}
	
	
	get.p <- function(null.row, obs.val){
		p <- length(which(null.row >= obs.val))/length(null.row)
		return(p)
		}
	
	all.p <- matrix(NA, ncol = 1, nrow = nrow(full.table))
	rownames(all.p) <- rownames(full.table)
	colnames(all.p) <- c("enrichment.p")
	for(i in 1:nrow(full.table)){
		all.p[i,1] <- get.p(null.dist[i,], full.table[i,1])
		}
	
	sig <- which(all.p[,1] <= sig.val)
	all.p[sig,,drop=FALSE]
	
	
	cols <- c("lightblue", "darkblue")
	col <- rep(cols, nrow(all.p))
	sig.idx <- (2*sig)-1
	col[sig.idx] <- "pink"
	
	# quartz(width = 10, height = 6)
	barplot(t(full.table), ylab = "Number of Genes", xlab = "Chromosome", beside = TRUE, col = col, main = plot.label)
	plot.dim <- par("usr")
	plot.height <- plot.dim[4] - plot.dim[3]
	par(xpd = TRUE)
	legend(x = 0, y = (plot.dim[4]+(plot.height*0.05)), fill = c(cols, "pink"), legend = c("Observed", "Expected", "Enriched"), horiz = TRUE)
	par(xpd = FALSE)

	chr.order <- order(suppressWarnings(as.numeric(query.table[,1])))
	sorted.table <- matrix(query.table[chr.order,], ncol = 1)
	rownames(sorted.table) <- rownames(query.table)

	invisible(list(sorted.table, full.table))
	}