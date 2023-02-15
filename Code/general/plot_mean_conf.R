#This function takes in a list of values and plots
#the mean for each element in the list surrounded 
#by a confidence interval estimate. It has an optional
#smoothing parameter set with spar.
#y.list is a list of values the same length as x
#each element in y corresponds to one eleent in x
#spar is a smoothing parameter between 0 and 1. 
#higher levels are more smooth

plot_mean_conf <- function(x, y.list, conf.level = 0.95, spar = 0,
  ylab = "", xlab = "", main = "", xlim = NULL, ylim = NULL){
  
  mean.y <- sapply(y.list, function(x) mean(x, na.rm = TRUE))
  conf.y <- t(sapply(y.list, function(x) t.test(x, conf.level = conf.level)$conf.int))

  smooth.mean <- smooth.spline(x, mean.y, spar = spar)
  smooth.upper <- smooth.spline(x, conf.y[,1], spar = spar)
  smooth.lower <- smooth.spline(x, conf.y[,2], spar = spar)

  plot.poly.xy(x, smooth.upper$y, x, smooth.lower$y, 
    border = NA, new.plot = TRUE, col = "gray",
    xlim = xlim, ylim = ylim)
  points(x, smooth.mean$y, col = "black", lwd = 3, type = "l")
  axis(1);axis(2)
  abline(v = c(0,1), lwd = 2, lty = 2, col = "gray")
  mtext(ylab, side = 2, line = 2.5)
  mtext(xlab, side = 1, line = 2.5)
  mtext(main, side = 3, line = 2.5)

  result <- rbind(smooth.upper$y, smooth.mean$y, smooth.lower$y)
  rownames(result) <- c("upper.conf", "mean", "lower.conf")
  colnames(result) <- x

  invisible(result)
}
