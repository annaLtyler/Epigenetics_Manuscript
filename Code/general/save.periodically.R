#This function is used for long processes that we want to 
#save periodically so that we don't lose too much progress.
#The function assumes that
save.periodically <- function(fun, max.iter, filename, save.every = 100){
    
    match.fun(fun)

    if(!file.exists(filename)){
        results <- vector(mode = "list", length = max.iter)
        for(i in 1:max.iter){
            results[[i]] <- fun(i)
            if(i %% save.every == 0){
                saveRDS(results, filename)
            }
        }
    }else{
        results <- readRDS(filename)
        next.pos <- which(sapply(results, length) == 0)
        
        #if there are no empty positions in the vector
        if(length(next.pos) == 0 && length(results) >= max.iter){
            cat("All iterations complete.\n")
            return()
        }else{
            if(length(next.pos) == 0 && length(results) < max.iter){
               next.pos <- length(results) + 1
            }else{
                next.pos <- min(next.pos)
            }
        }

        for(i in next.pos:max.iter){
            results[[i]] <- fun(i)
            if(i %% save.every == 0){
                saveRDS(results, filename)
            }
        }
    }
}