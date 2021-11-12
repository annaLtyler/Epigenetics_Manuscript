#This function centers data based on given batches
#such that all the batches have the same mean value.

batch.center <- function(gene.expr, batch){

    batch.levels <- unique(batch)
    batch.idx <- sapply(batch.levels, function(x) which(batch == x))
    batch.means <- sapply(batch.idx, function(x) mean(gene.expr[x], na.rm = TRUE))
    not.na.locale <- which(!is.na(batch.means))
    if(length(not.na.locale) < 2){
        return(gene.expr)
    }else{
        new.mean <- mean(batch.means, na.rm = TRUE)
        mean.diff <- batch.means - new.mean
        adj.expr <- gene.expr
        for(i in not.na.locale){
            adj.expr[batch.idx[[i]]] <- gene.expr[batch.idx[[i]]] - mean.diff[i]
        }
    }

    return(adj.expr)
}
