#This function is for shifting columns of the methylation
#table to see if we can remove the cascading effect

adjust.methyl.table <- function(methyl.table, col.num, num.spaces){
	
	new.table <- methyl.table

	for(i in 1:length(col.num)){
		padding <- rep(0, abs(num.spaces[i]))
		
		if(num.spaces[i] > 0){
			new.col <- c(padding, methyl.table[,col.num[i]])
			}
		if(num.spaces[i] < 0){
			new.col <- c(methyl.table[length(padding):nrow(methyl.table),col.num[i]], padding)
			}
		
		new.table[,col.num[i]] <- new.col[1:nrow(methyl.table)]
		# cbind(methyl.table[,col.num[i]], new.table[,col.num[i]])
		}

	return(new.table)
}