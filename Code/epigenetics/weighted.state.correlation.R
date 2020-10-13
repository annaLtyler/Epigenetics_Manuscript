#This function looks at whether combinations of weighted states 
#improve the correlation between chromatin state proportion and 
#both inbred expression and DO cis eQTL coefficient

weighted.state.correlation <- function(DO.results.file, inbred.results.file, 
do.cis.coef, chrom.state.props, col.table, strain.mean.expr, state.weights){

    #match up all the data objects
    matched.data <- sort.data(do.cis.coef, chrom.state.props, strain.mean.expr, 
    col.table)
    common.prop <- matched.data$chromatin.prop
    common.coef <- matched.data$cis.DO
    common.expr <- matched.data$inbred.expr

    #scale all expression values
    scaled.expr <- lapply(common.expr, scale)
    
    state.pairs <- pair.matrix(1:length(state.weights), self.pairs = TRUE)

    if(!file.exists(DO.results.file) || !file.exists(inbred.results.file)){
        base.weights <- rep(1, ncol(common.prop[[1]]))
        all.DO.cor <- all.inbred.cor <- list()
        
        for(i in 1:nrow(state.pairs)){
        
            one.state <- rep(0, length(base.weights))
            weight.idx <- unique(c(state.pairs[i,1], state.pairs[i,2]))
            one.state[weight.idx] <- state.weights[weight.idx]
            weighted.states <- lapply(common.prop, function(x) apply(x, 1, function(y) y*one.state))
            weighted.sums <- lapply(weighted.states, function(x) colSums(x, na.rm = TRUE))
            
            strain.locale <- order.strains(names(weighted.sums[[1]]), 
            names(common.coef[[1]]), col.table)

            #correlate weighted states with DO coefficients
            all.DO.cor[[i]] <- mapply(function(x,y) 
            cor(x, as.vector(y)[strain.locale], use = "complete"), 
            weighted.sums, common.coef)

            #correlate weighted states with scaled expression
            all.inbred.cor[[i]] <- mapply(function(x,y) 
            cor(x, as.vector(y), use = "complete"), weighted.sums, scaled.expr)
        }

        #for the final test use all weights
        weighted.states <- lapply(common.prop, function(x) apply(x, 1, function(y) y*state.weights))
        weighted.sums <- lapply(weighted.states, function(x) colSums(x, na.rm = TRUE))
        
        all.DO.cor[[(i+1)]] <- mapply(function(x,y) 
        cor(x, as.vector(y)[strain.locale], use = "complete"), weighted.sums, common.coef)

        all.inbred.cor[[(i+1)]] <- mapply(function(x,y) 
        cor(x, as.vector(y), use = "complete"), weighted.sums, scaled.expr)

        names(all.DO.cor) <- names(all.inbred.cor) <- c(apply(state.pairs, 1, function(x) paste(unique(x), collapse = "_")), "all")
        saveRDS(all.DO.cor, DO.results.file)
        saveRDS(all.inbred.cor, inbred.results.file)

    }else{
        all.DO.cor <- readRDS(DO.results.file)
        all.inbred.cor <- readRDS(inbred.results.file)
    }

    return(list("DO.coef" = all.DO.cor, "inbred.expr" = all.inbred.cor))
        
}