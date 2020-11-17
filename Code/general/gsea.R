#Performs GSEA as described in Subramanian et al. 2005
#ordered.genes is an ordered numeric vector with gene
#name labels. go.genes is a vector of genes in the current
#GO term (from read.gmt) or another gene list you are 
#interested in. min.size and max.size are the minimum 
#and maximum gene list size to consider. nperm is the
#number of permutations to perform. p is the weighting 
#factor. If p = 0, the enrichment score reduces to the 
#Kolmogorov-Smirnov test stastic. When p = 1, genes are
#weighted by their numeric value.

gsea <- function(ordered.genes, go.genes, min.size = 15, max.size = 600,
nperm = 100, p = 0){
    
    n.total.genes <- length(ordered.genes)
    n.set.genes <- length(go.genes)
    n.intersect <- length(intersect(names(ordered.genes), go.genes))

    if(n.intersect < min.size || n.set.genes > max.size){
        result <- list("Enrichment" = NA)
        return(result)
        }
    
    get_enrich <- function(hit.idx){
        hits.to.use <- sapply(1:n.total.genes, function(x) hit.locale[which(hit.idx < x)])
        weighted.hits <- sapply(hits.to.use, function(x) sum(abs(ordered.genes[x])^p, na.rm = TRUE))
        phit <- weighted.hits/weighted.hits[n.total.genes]
        #plot(phit)
    
        n.hits <- sapply(hits.to.use, function(x) length(x))
        weighted.misses <- (1:n.total.genes) - n.hits
        pmiss <- weighted.misses / weighted.misses[n.total.genes]
        #plot(pmiss)

        enrich = phit - pmiss
        #plot(enrich, type = "l");abline(h = 0)
        return(enrich)
    }

    #where true hits are
    hit.locale <- which(names(ordered.genes) %in% go.genes)
    true.enrich <- get_enrich(hit.locale)
    #plot(true.enrich, type = "l");abline(h = 0)
    true.ES <- max(true.enrich)

    #permute location of hits
    if(nperm > 0){
        perm.ES <- rep(NA, nperm)
        for(p in 1:nperm){
            perm.locale <- sample(1:length(ordered.genes), length(hit.locale))    
            perm.enrich <- get_enrich(perm.locale)
            #plot(perm.enrich, type = "l");abline(h = 0)
            perm.ES[p] <- max(perm.enrich)
        }
        #hist(perm.ES);abline(v = true.ES, col = "red")

    result <- list("Enrichment" = true.enrich, "Perm.ES" = perm.ES)
    return(result)
    }else{
        result <- list("Enrichment" = true.enrich)
    }
    
    return(result)   
}