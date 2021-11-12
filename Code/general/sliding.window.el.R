#This function returns a lits of elements to use
#for a sliding window

sliding.window.el <- function(elV, window.size, gap.size, plot.results = FALSE){
	
	total.num = length(elV)
	el.list <- list()
	
	start.pos <- 1
	list.ind <- 1
	while(start.pos + window.size <= total.num){
		el.list[[list.ind]] <- elV[start.pos:(start.pos+window.size-1)]
		list.ind = list.ind + 1
		start.pos = start.pos + gap.size
		}
	if(start.pos < total.num){
		el.list[[list.ind]] <- elV[start.pos:total.num]
		}

	if(plot.results){
		plot.new()
		plot.window(xlim = c(min(elV), max(elV)), ylim = c(1, length(el.list)))
		for(i in 1:length(el.list)){
			points(x = el.list[[i]], y = rep(i, length(el.list[[i]])))
		}
		axis(1)
	}


	return(el.list)
	
	
	
	
}
