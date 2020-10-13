#This function gets the mean methylation values
#around the center point of values that have been
#centered by center.on.feature().
#The upstream and downstream buffer are in bp
#if is.region is set to TRUE, the function instead
#finds three means: the mean for positions below 0,
#for positions between 0 and 1, and for positions
#above 1.

get.mean.near.center <- function(vals, is.region = FALSE, upstream.buffer = 5000, 
downstream.buffer = 5000, get.mean = TRUE){
  
  if(length(vals) == 1 && is.na(vals)){
    if(is.region){
      return(rep(NA, 3))
    }else{
      return(NA)
    }
  }

  pos <- as.numeric(names(vals))

  if(is.region){
    up.locale <- which(pos < 0)
    mid.locale <- intersect(which(pos >= 0), which(pos <= 1))
    down.locale <- which(pos > 1)
    up.mean <- mean(vals[up.locale], na.rm = TRUE)
    mid.mean <- mean(vals[mid.locale], na.rm = TRUE)
    down.mean <- mean(vals[down.locale], na.rm = TRUE)
    result <- c("up" = up.mean, "mid" = mid.mean, "down" = down.mean)
    return(result)
  }else{

    center.locale <- get.nearest.pt(pos, 0)
    center.pos <- pos[center.locale]

    #plot(as.numeric(names(vals)), vals)
    #abline(v = c(center.pos, -upstream.buffer, downstream.buffer))
    

    if(abs(center.pos) > upstream.buffer && abs(center.pos) > downstream.buffer){
      return(NA)
    } 

    min.pos <- 0 - upstream.buffer
    max.pos <- 0 + downstream.buffer
    within.pos <- intersect(which(pos >= min.pos), which(pos <= max.pos))
    if(length(within.pos) > 0){
      if(get.mean){
        return(mean(vals[within.pos], na.rm = TRUE))
      }else{
        return(vals[within.pos])
      }
    }else{
      return(NA)
    }
  }
}
