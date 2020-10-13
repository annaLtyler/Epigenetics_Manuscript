#This function turns a SNP matrix into a 
#tree
snps2tree <- function(snp.mat){
	require(ape)
	require(e1071)
	
	u_snps <- unique(unlist(apply(snp.mat, 2, unique)))
	key <- 1:length(u_snps)
	
	keyed.mat <- snp.mat
	for(i in 1:length(u_snps)){
		keyed.mat[which(snp.mat == u_snps[i])] <- key[i]
		}
	
	
	keyed.mat <- apply(keyed.mat, 2, as.numeric)
	
	dist.mat <- matrix(0, nrow = ncol(snp.mat), ncol = ncol(snp.mat))
	col.pairs <- pair.matrix(1:ncol(snp.mat))
	for(i in 1:nrow(col.pairs)){
		h.dist <- hamming.distance(keyed.mat[,col.pairs[i,1]],keyed.mat[,col.pairs[i,2]])
		dist.mat[col.pairs[i,1], col.pairs[i,2]] <- h.dist
		dist.mat[col.pairs[i,2], col.pairs[i,1]] <- h.dist
		}
	
	colnames(dist.mat) <- rownames(dist.mat) <- colnames(snp.mat)
	dist.tree <- nj(dist.mat)
	# plot(dist.tree, main = main)
	return(dist.tree)
}