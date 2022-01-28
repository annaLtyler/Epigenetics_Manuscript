#This function gets lists of genes with the specified state 
#at binned positions along the gene body and for each strain. 
#The enrichments are gene lists are saved as filename, and
#only calculated if the file does not exist.
#The enrichments are then processed to compare enrichments
#across strains.
#to plot this, you need to geneate a very wide device, like
#20 inches.

state_peak_genes <- function(binned.chromatin.states, state.id, min.n = 30, enrich.filename, 
    gene.filename, max.term.size = NULL, n.terms = 10, num.strain = 9, 
    sort.by = c("default", "p_value")){
    

    if(!file.exists(enrich.filename)){
        one.state <- lapply(binned.chromatin.states, function(x) if(length(x) > 1){t(x[[state.id]])})
        not.null.chrom <- which(sapply(one.state, length) > 0)

        #return 
        state_in_gene <- function(state.chrom.mat){    
            strains.which <- apply(state.chrom.mat, 1, function(x) length(which(x == 1)))
            pos.which <- apply(state.chrom.mat, 2, function(x) length(which(x == 1)))
            if(length(pos.which) == 0){
                return(NULL)
            }else{
                pos.with.state <- state.chrom.mat
                pos.with.state[pos.which]  <- 1
                return(pos.with.state)
            }
        }

        pos_enrich <- function(strain.list, min.n){
            list.size <- sapply(strain.list, length)
            big.enough <- which(list.size >= min.n)
            enrich.list <- lapply(strain.list[big.enough], function(x) gost(names(x), organism = "mmusculus", sources = c("GO", "KEGG", "REACTOME")))
            names(enrich.list) <- round(as.numeric(names(strain.list)[big.enough]), 2)
            #plot.enrichment.group(enrich.list)
            return(enrich.list)
        }

        #split the chromatin state matrices by strain
        state.by.strain <- lapply(1:num.strain, function(x) sapply(not.null.chrom, function(y) one.state[[y]][x,]))
        names(state.by.strain)  <- rownames(one.state[[not.null.chrom[1]]])

        #for each position in each matrix, get the list of genes present in that position
        state.by.pos <- lapply(state.by.strain, function(x) apply(x, 1, function(y) which(y == 1)))
        saveRDS(state.by.pos, gene.filename)
        #barplot(sapply(state.by.pos[[1]], length))
        
        #get enrichments for each set of genes greater than the minimum set in the arguments
        
        enrich.by.strain.pos <- lapply(state.by.pos, function(x) pos_enrich(x, min.n))
        saveRDS(enrich.by.strain.pos, enrich.filename)
    }else{
        state.by.pos <- readRDS(gene.filename)
        enrich.by.strain.pos <- readRDS(enrich.filename)
    }

    #count genes per strain with mark at each position
    gene.count.mat <- sapply(state.by.pos, function(x) sapply(x, length))
    #pheatmap(gene.count.mat, cluster_rows = FALSE, cluster_cols = FALSE, scale = "column")

    #compare enrichments by position across strains
    #make enrichment matrices for each strain
    enrich.mat.by.strain <- lapply(enrich.by.strain.pos, 
        function(x) plot.enrichment.group(x, plot.results = FALSE, max.term.size = max.term.size,
        n.terms = n.terms, sort.by = sort.by))
    #plot.enrichment.group(enrich.by.strain.pos[[4]], cluster_cols = FALSE)

    unique.terms <- unique(unlist(lapply(enrich.mat.by.strain, rownames)))
    unique.pos <- unique(unlist(lapply(enrich.mat.by.strain, colnames)))

    #Make an array of common matrices to look across all strains
    enrich.array <- array(NA, dim = c(length(unique.terms), length(unique.pos), num.strain))
    dimnames(enrich.array) <- list(unique.terms, sort(as.numeric(unique.pos)), names(state.by.pos))

    for(strain in 1:length(enrich.mat.by.strain)){
        combined.mat <- count.mat <- matrix(0, nrow = length(unique.terms), ncol = length(unique.pos))
        rownames(combined.mat) <- unique.terms
        colnames(combined.mat) <- sort(as.numeric(unique.pos))

        sig.terms <- which(enrich.mat.by.strain[[strain]] > 0, arr.ind = TRUE)
        if(length(sig.terms) > 0){
            for(j in 1:nrow(sig.terms)){
                term.name <- rownames(enrich.mat.by.strain[[strain]])[sig.terms[j,1]]
                pos.name <- colnames(enrich.mat.by.strain[[strain]])[sig.terms[j,2]]
                combined.mat[term.name, pos.name] <- enrich.mat.by.strain[[strain]][sig.terms[j,1], sig.terms[j,2]]
            }
            enrich.array[,,strain] <- combined.mat
        }
    }

    #pheatmap(enrich.array[,,4], cluster_cols = FALSE) 
    all.pos.vals <- t(sapply(1:dim(enrich.array)[3], function(x) colSums(enrich.array[,,x])))
    all.pos.sum <- colSums(all.pos.vals)
    #collect terms seen at each position

    extract_terms <- function(array.mat){
        term.log.p <- rowSums(array.mat[which(rowSums(array.mat) > 0),,drop=FALSE])
        if(length(term.log.p) == 0){
        term.log.p <- 0
        names(term.log.p) <- "no enrichment"
        }
        return(term.log.p)
    }
    
    term.by.pos <- lapply(1:dim(enrich.array)[2], function(x) extract_terms(enrich.array[,x,]))
    #word.bags <- lapply(term.by.pos, function(x) sort(table(str_to_lower(unlist(strsplit(x, " ")))), decreasing = TRUE))
    term.bags <- lapply(term.by.pos, function(x) sort(unlist(x), decreasing = TRUE))

    plot.file <- paste0("state.position.enrichment.", state.id, ".pdf")
    #pdf(file.path("~/Desktop", plot.file), height = 6, width = 20)
    par(mar = c(4,4,2,70))
    a <- barplot(all.pos.vals, col = col.table[,3], las = 2, horiz = TRUE, axes = FALSE, cex.names = 0.7)
    word.y <- segment.region(min(a), max(a), length(term.bags), "ends")    

    par(xpd = NA)
    for(i in 1:length(term.bags)){
        text(x = all.pos.sum[i], y = word.y[i], labels = paste(names(term.bags[[i]]), collapse = ", "), cex = 0.5, adj = 0)
    }
    par(xpd = FALSE)
    mtext(paste("State", state.id), side = 3, outer = TRUE, line = -2.5)
    #dev.off()

    num.pos <- length(all.pos.sum)
    genes.by.pos <- lapply(1:num.pos, function(x) unique(unlist(lapply(state.by.pos, function(y) names(y[[x]])))))
    names(genes.by.pos) <- names(all.pos.sum)
    #collapse the enrichment terms across strains.
    collapsed.array <- Reduce("+", lapply(1:dim(enrich.array)[3], function(x) enrich.array[,,x]))

    result <- list("genes.by.position" = genes.by.pos, 
    "enrichment.by.position" = collapsed.array, "gene.counts" = gene.count.mat)
    invisible(result)


}