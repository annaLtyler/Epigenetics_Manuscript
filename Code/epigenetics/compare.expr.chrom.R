#This function calculates whether a chromatin distance matrix 
#is correlated with inbred expression
	
compare.expr.chrom <- function(d.mat, expr.list){
	non.zero <- which(d.mat != 0)
	na.locale <- which(is.na(d.mat))
	if(length(non.zero) > 0 && length(na.locale) == 0){
		fit <- cmdscale(d.mat, k = 1)
		test <- aov(unlist(lapply(expr.list, median))~fit)
		expr.state.pval <- anova(test)$"Pr(>F)"[1]
		return(expr.state.pval)
		}else{
		return(NA)	
		}
	}
