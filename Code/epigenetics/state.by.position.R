#This function calculates state presence/absence across the set
#of given genes. The gene.ids need to be of the same form as the
#names of binned.chrom. binned.chrom is a list
#of binned chromatin matrices from bin.centered.chromatin().
#This function calculates the amount of the given state across
#all genes as a function of the binned chromatin positions. 

state.by.position <- function(gene.ids, group.gene.expr, state.id, binned.chrom){

    #get the binned chromatin matrices for these genes as well.
    gene.locale <- match(gene.ids, names(binned.chrom))
    gene.chrom.bins <- binned.chrom[gene.locale]

    one.state <- lapply(gene.chrom.bins, function(x) if(length(x) > 1){t(x[[state.id]])})
    not.null.chrom <- which(sapply(one.state, length) > 0)
    not.null.expr <- which(sapply(group.gene.expr, length) > 1)
    not.null <- intersect(not.null.chrom, not.null.expr)
    
    state.sum.mat <- Reduce("add.mat", one.state[not.null])
    #pheatmap(state.sum.mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE)
    #plot.new()
    #plot.window(xlim = c(-2, 2), ylim = c(0, max(state.sum.mat, na.rm = TRUE)))
    #for(i in 1:nrow(state.sum.mat)){
    #    points(as.numeric(colnames(state.sum.mat)), state.sum.mat[i,], type = "l", col = col.table[i,3], lwd = 3)
    #}
    #axis(1);abline(v = c(0,1))

    state.count.mat <- Reduce("+", lapply(one.state[not.null], state.to.count))

    #plot.new()
    #plot.window(xlim = c(-2, 2), ylim = c(0, max(state.count.mat, na.rm = TRUE)))
    #for(i in 1:nrow(state.sum.mat)){
    #    points(as.numeric(colnames(state.count.mat)), state.count.mat[i,], type = "l", col = col.table[i,3], lwd = 3)
    #}
    #axis(1);abline(v = c(0,1))

    state.ssq.mat <- Reduce("ssq.mat", one.state[not.null])

    state.sum.vec <- colSums(state.sum.mat, na.rm = TRUE)
    state.count.vec <- colSums(state.count.mat, na.rm = TRUE)
    state.avg <- state.sum.vec/state.count.vec
    #plot(state.avg, type = "l")
    state.var <- colSums(state.ssq.mat, na.rm = TRUE)/state.count.vec - (state.avg)^2
    state.sd <- sqrt(state.var)
    state.se <- state.sd/sqrt(state.count.vec)

    results <- list("position" = as.numeric(colnames(state.sum.mat)), 
        "state.avg" = state.avg, "state.sd" = state.sd, "state.se" = state.se)
    return(results)
}
