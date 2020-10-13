#returns difference in mean expression between
#groups. The groups can be a list with multiple
#sub groups, for example, strains in control v treatment
diff.expr <- function(rna.seq, idx, group1, group2){
	
	if(class(group1) == "list"){
		group1.expr <- lapply(group1, function(x) rna.seq[idx, x])
		group2.expr <- lapply(group2, function(x) rna.seq[idx, x])
		group.diffs <- vector(mode = "list", length = length(group1))
		names(group.diffs) <- names(group1)
		for(i in 1:length(group.diffs)){
			group.diffs[i] <- mean(group1.expr[[i]], na.rm = TRUE) - mean(group2.expr[[i]], na.rm = TRUE)
			}
		result <- unlist(group.diffs)
		}else{
		result <- mean(rna.seq[idx, group1], na.rm = TRUE) - mean(rna.seq[idx,group2], na.rm = TRUE)
		}
	
	return(result)
	
}