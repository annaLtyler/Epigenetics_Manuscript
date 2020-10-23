#This function calculates state presence/absence across the set
#of given genes. The gene.ids need to be of the same form as the
#names of binned.chrom. binned.chrom is a list
#of binned chromatin matrices from bin.centered.chromatin().
#This function calculates the amount of the given state across
#all genes as a function of the binned chromatin positions. 

state.by.position <- function(gene.ids, state.id, binned.chrom){

    #get the binned chromatin matrices for these genes as well.
    gene.locale <- match(gene.ids, names(binned.chrom))
    gene.chrom.bins <- binned.chrom[gene.locale]

    one.state <- lapply(gene.chrom.bins, function(x) if(length(x) > 1){t(x[[state.id]])})
    not.null.chrom <- which(sapply(one.state, length) > 0)
    not.null.expr <- which(sapply(group.gene.expr, length) > 1)
    not.null <- intersect(not.null.chrom, not.null.expr)
    #one.state.mat <- Reduce("rbind", one.state[not.null])
    #pheatmap(one.state.mat, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)
    
    state.sum.mat <- Reduce("add.mat", one.state[not.null])
    state.count.mat <- Reduce("+", lapply(one.state[not.null], state.to.count))
    state.ssq.mat <- Reduce("ssq.mat", one.state[not.null])

    state.sum.vec <- colSums(state.sum.mat, na.rm = TRUE)
    state.count.vec <- colSums(state.count.mat, na.rm = TRUE)
    state.avg <- state.sum.vec/state.count.vec
    state.var <- colSums(state.ssq.mat, na.rm = TRUE)/state.count.vec - (state.avg)^2
    state.sd <- sqrt(state.var)
    state.se <- state.sd/sqrt(state.count.vec)

    results <- list("position" = as.numeric(colnames(state.sum.mat)), 
        "state.avg" = state.avg, "state.sd" = state.sd, "state.se" = state.se)
    return(results)
}
