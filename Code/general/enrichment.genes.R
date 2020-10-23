#This function extracts gene names associated with
#specific terms from a gprofiler enrichment list.
#if this is a group of enrichment tables, set 
#is.group to TRUE. If it is just enrichments for 
#a single list, set is.group to FALSE.
#to run this function, evcodes must be set to TRUE
#when running gost.

enrichment.genes <- function(enrichment.list, term.name, is.group = TRUE){

    if(is.group){
        enrichment.tables <- lapply(enrichment.list, function(x) x$result)
    }else{
        enrichment.tables <- list(enrichment.list$result)
    }

    intersection.column.exists <- which(colnames(enrichment.tables[[1]]) == "intersection")
    if(length(intersection.column.exists) == 0){
        stop("evcodes must be set to TRUE when running gost.")
    }

    all.terms <- lapply(enrichment.tables, function(x) x[,"term_name"])
    term.locale <- grep(term.name, all.terms)
    if(length(term.locale) == 0){
        cat("I can't find the term. Check spelling and capitalization.\n")
        return(NULL)
    }
    all.genes <- lapply(enrichment.tables, function(x) x[,"intersection"])

    gene.list <- vector(mode = "list", length = length(term.locale))
    names(gene.list) <- names(enrichment.list)[term.locale]

    for(i in 1:length(term.locale)){
        term.idx <- grep(term.name, all.terms[[term.locale[i]]])
        gene.list[[i]] <- strsplit(all.genes[[term.locale[i]]][term.idx], ",")[[1]]
    }

    return(gene.list)
}