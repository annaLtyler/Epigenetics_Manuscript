#This function uses kmeans clustering to cluster
#a matrix, either as is, or using its svd with
#a given number of pcs.
#It returns a membership vector indicating the
#cluster assignment for each row in the matrix.

kmeans_cluster_mat <- function(mat, min.cl = 2, max.cl = 20, use.pc = TRUE, pc = 2, 
    plot.results = FALSE){


    decomp <- plot.decomp(mat, pc = pc, plot.results = FALSE)

    if(use.pc){
        k.test <- test.pam.k(decomp$u, min.cl:max.cl, plot.results = FALSE)
    }else{
        k.test <- test.pam.k(mat, min.cl:max.cl, plot.results = FALSE)
    }

    mean.cl <- sapply(k.test[[1]], mean)
    best_clustering <- which.max(mean.cl)
    best_cl <- k.test[[2]][,best_clustering]
    u_cl <- unique(best_cl)    
    #cluster membership
    cl_mem <- lapply(u_cl, function(x) which(best_cl == x))
        
    #plot PC with clusters colored
    if(plot.results){
        pairs(decomp$u, col = best_cl, pch = 16)
    }

    return(best_cl)

}
