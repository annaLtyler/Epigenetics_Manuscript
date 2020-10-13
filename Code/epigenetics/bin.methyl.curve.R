#This is a special version of bin.curve that
#groups methylation curves into islands.
#Instead of binning peak to peak, it looks
#for sharp jumps in methylation percent
#and bins at the changes rather than peak to
#trough. transition.percent tells the function
#where island blocks switch from being methylated
#to unmethylated.

bin.methyl.curve <- function(the.curve, plot.peaks = FALSE, transition.percent = 50){
	
    if(plot.peaks){
        plot(the.curve, type = "l")
        abline(h = transition.percent)
    }

    #smooth out the curve using transition.percent
    flattened.curve <- the.curve
    flattened.curve[which(the.curve <= transition.percent)] <- 0
    flattened.curve[which(the.curve > transition.percent)] <- 100
    if(plot.peaks){
        points(flattened.curve, type = "l", col = "blue")
    }

    consec.pts <- consec.pairs(flattened.curve)
    transition.side <- consec.pts > transition.percent
    transition.pts <- which(apply(transition.side, 1, function(x) x[1] != x[2])) + 1

    if(plot.peaks){
        points(y = rep(transition.percent, length(transition.pts)), 
        x = transition.pts+0.5, pch = "*", cex = 1.5, col = "red")
    }

#    consec.transitions <- consec.pairs(transition.pts)
#    block.size <- consec.transitions[,2] - consec.transitions[,1]
#    small.blocks <- which(block.size < min.block.size)

#    while(length(small.blocks) > 0){
#        remove.transition <- small.blocks[1]+1
#        transition.pts <- transition.pts[-remove.transition]
#        consec.transitions <- consec.pairs(transition.pts)
#        block.size <- consec.transitions[,2] - consec.transitions[,1]
#        small.blocks <- which(block.size < min.block.size)
#    }
    
    bins <- rep(1,length(the.curve))
    bin.num <- 2
    for(i in 1:(length(transition.pts)-1)){
        start.bin <- transition.pts[i]
        end.bin <- transition.pts[i+1]
        bins[start.bin:end.bin] <- bin.num
        bin.num <- bin.num + 1
    }
    bins[end.bin:length(bins)] <- bin.num

    return(bins)
}