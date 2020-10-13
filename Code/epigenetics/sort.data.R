#This function takes the DO cis eQTL object, the chromatin state proportions object,
#and strain mean expression. It sorts by gene name and by strain and returns three 
#sorted lists


sort.data <- function(do.cis.coef, chrom.state.props, strain.mean.expr, col.table){

    #put the chrom.prop matrices in the same strain order as expression.
    do.strains <- names(all.cis.coef[[1]])
    chrom.strain.order <- order.strains(do.strains, colnames(chrom.state.props[[1]]), col.table)
    chrom.strain.order <- chrom.strain.order[which(!is.na(chrom.strain.order))]
    ordered.props <- lapply(chrom.state.props, function(x) 
    if(length(x) > 1){t(x[,chrom.strain.order])}else{NA})

    expr.order <- order.strains(do.strains, names(strain.mean.expr[[1]]), col.table)
    expr.order <- expr.order[which(!is.na(expr.order))]
    #order expression values
    ordered.expr <- lapply(strain.mean.expr, function(x) x[expr.order])
    
    #make sure all data objects are listing the same genes in the same order
    common.genes <- Reduce("intersect", list(names(ordered.props), names(ordered.expr), 
    names(all.cis.coef)))
    common.prop.locale <- match(common.genes, names(ordered.props))
    common.coef.locale <- match(common.genes, names(all.cis.coef))
    common.expr.locale <- match(common.genes, names(strain.mean.expr))

    common.prop <- ordered.props[common.prop.locale]
    common.coef <- all.cis.coef[common.coef.locale]
    common.expr <- ordered.expr[common.expr.locale]

    #some of the proportion matrices are not available. 
    #Remove these indices from all lists.
    has.prop.vals <- which(sapply(common.prop, length) > 1)
    has.expr.vals <- which(sapply(common.expr, function(x) !is.na(x[1])))
    has.vals <- intersect(has.prop.vals, has.expr.vals)

    common.prop <- common.prop[has.vals]
    common.coef <- common.coef[has.vals]
    common.expr <- common.expr[has.vals]


    results <- list("chromatin.prop" = common.prop, "cis.DO" = common.coef, 
    "inbred.expr" = common.expr)
    return(results)

}