#This function calculates state presence/absence across the set
#of given genes. The gene.ids need to be of the same form as the
#names of binned.chrom and scaled.expr. binned.chrom is a list
#of binned chromatin matrices from bin.centered.chromatin().
#scaled.expr is a list of scaled expression values by gene.
#This function calculates the amount of the given state across
#all genes as a function of the binned chromatin positions. 
#It also calculates the correlation between gene expression 
#across strains and the chromatin state in each position.

state.expr.by.position <- function(gene.ids, state.id, binned.chrom, scaled.expr){
    #get the scaled expression for the gene ids
    gene.locale <- match(gene.ids, names(scaled.expr))
    group.gene.expr <- scaled.expr[gene.locale]
    exprV <- unlist(lapply(group.gene.expr, function(x) x[,1]))

    #get the binned chromatin matrices for these genes as well.
    gene.locale <- match(gene.ids, names(binned.chrom))
    gene.chrom.bins <- binned.chrom[gene.locale]

    one.state <- lapply(gene.chrom.bins, function(x) if(length(x) > 1){t(x[[state.id]])})
    not.null.chrom <- which(sapply(one.state, length) > 0)
    not.null.expr <- which(sapply(group.gene.expr, length) > 1)
    not.null <- intersect(not.null.chrom, not.null.expr)
    #one.state.mat <- Reduce("rbind", one.state[not.null])
    #pheatmap(one.state.mat, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)
    pheatmap(one.state.mat[order(exprV),], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)
    
    diff.state.mat <- one.state.mat
    diff.state.mat[which(is.na(one.state.mat))] <- 0
    diff.state.mat[which(one.state.mat == 1)] <- 1
    diff.state.mat[which(one.state.mat == 0)] <- -1
    
    pheatmap(diff.state.mat, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)
    one.state.step <- diff.state.mat[,1,drop=FALSE]
    for(i in 2:ncol(diff.state.mat)){
        one.state.step <- cbind(one.state.step, (one.state.step[,(i-1)] + diff.state.mat[,i]))
    }
    pheatmap(one.state.step, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)
    pheatmap(one.state.step[order(exprV),], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)
    
    state.sum.mat <- Reduce("add.mat", one.state[not.null])
    state.count.mat <- Reduce("+", lapply(one.state[not.null], state.to.count))
    state.ssq.mat <- Reduce("ssq.mat", one.state[not.null])

    state.sum.vec <- colSums(state.sum.mat, na.rm = TRUE)
    state.count.vec <- colSums(state.count.mat, na.rm = TRUE)
    state.avg <- state.sum.vec/state.count.vec
    state.var <- colSums(state.ssq.mat, na.rm = TRUE)/state.count.vec - (state.avg)^2
    state.se <- sqrt(state.var)/sqrt(state.count.vec)

    #multiply state matrices by expression
    state.expr <- lapply(not.null, function(x) diag(group.gene.expr[[x]][,1]) %*% one.state[[x]])
    state.expr.sum <- Reduce("add.mat", state.expr)
    state.expr.avg <- colSums(state.expr.sum, na.rm = TRUE)/state.count.vec
    state.expr.cor <- state.expr.avg/sqrt(state.var)

    results <- list("position" = as.numeric(colnames(state.sum.mat)), 
        "state.avg" = state.avg, "state.se" = state.se, "expr.avg" = state.expr.avg,
        "expr.cor" = state.expr.cor)
    return(results)
}
