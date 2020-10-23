#This function correlates a gene's binned state matrix with expression
#cent.chrom is the chromatin matrix with coordinates centered on the 
#TSS and TES, gene.expr is the scaled expression of the gene.

binned.state.expr.cor <- function(cent.chrom, gene.expr, col.table, bin.min = 0, bin.max = 1,
num.states = 8){

    results <- list("r" = rep(NA, num.states), "p" = rep(NA, num.states))
    
    if(all(is.na(gene.expr))){
        return(results)
    }

    if(length(unique(as.vector(cent.chrom))) < 2){
        return(results)
    }

    pos <- as.numeric(rownames(cent.chrom))
    bin.coord <- intersect(which(pos >= bin.min), which(pos <= bin.max))

    if(length(bin.coord) == 0){
        return(results)
    }


    get.state.prop <- function(v, s){
        state.prop <- length(which(v == s))/length(v)
        return(state.prop)
    }

    binned.chrom <- cent.chrom[bin.coord,,drop=FALSE]
    state.strain.props <- apply(binned.chrom, 2, function(y) sapply(1:num.states, function(b) get.state.prop(y, b)))
    rownames(state.strain.props) <- 1:num.states

    strain.order <- match.order(colnames(state.strain.props), rownames(gene.expr), col.table)
    cor.tests <- suppressWarnings(apply(state.strain.props, 1, function(x) cor.test(gene.expr[strain.order,1], x)))
    cor.r <- sapply(cor.tests, function(x) x$estimate)
    cor.p <- sapply(cor.tests, function(x) x$p.value)

    results <- list("r" = cor.r, "p" = cor.p)
    return(results)
}