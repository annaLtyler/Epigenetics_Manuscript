#This function is for keeping strain orders 
#consistent and a bit simpler.
#given two vectors of strain names, it returns
#the order the second vector should be put in
#to match the first. It expects col.table, which
#is stored in the file strain.color.table.txt
#If this file is changed, this function will need
#to change

order.strains <- function(strain.fixed, strain.reorder, col.table){
	
	#figure out which column each strain vector matches the best
	match.col <- function(strainV, col.table){
		strain.matches <- apply(col.table, 2, function(x) match(strainV, x))
		best.match <- which.max(apply(strain.matches, 2, function(x) length(which(!is.na(x)))))
		if(length(best.match) == 0){stop("Couldn't find a match for the strain names:", strainV, "\n")}
		return(as.vector(best.match))
		}
	
	strain.fixed.col <- match.col(strain.fixed, col.table)
	strain.reorder.col <- match.col(strain.reorder, col.table)
	
	#sort col.table to match the order of the strain.fixed vector
	fixed.order <- match(strain.fixed, col.table[,strain.fixed.col])
	#cbind(strain.reorder, col.table[fixed.order,])
	fixed.order.table <- col.table[fixed.order,]

	#Then reorder the strain.reorder vector to match that vector
	reorder.order <- match(fixed.order.table[,strain.reorder.col], strain.reorder)
	#cbind(strain.reorder[reorder.order], strain.fixed)

	return(reorder.order)

}