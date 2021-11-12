#This function plots values that have been centered around
#a feature by center.on.feature() and then processed by 
#plot.centered.vals. Instead of a list, this function takes
#in the matrix returned by plot.centered.vals. It is 
#a shortcut to just plotting the values without the 
#need to align them. 

plot.centered.mat <- function(val.mat, min.representation = 10, 
plot.fun = c("mean", "median"), plot.label = "Center", show.var = c("se", "sd"),
plot.individual = FALSE, ylim = c(0, 50), ylab = "Methylation Percent Mean", 
plot.line = TRUE, plot.hex = FALSE, min.upstream = -2, max.downstream = 2){
  
  plot.fun <- plot.fun[1]
  show.var = show.var[1]

  all.pos <- as.numeric(colnames(val.mat))
  min.x <- max(c(floor(min(all.pos, na.rm = TRUE)), min.upstream))
  max.x <- min(c(ceiling(max(all.pos, na.rm = TRUE)), max.downstream))

  min.val <- min(val.mat, na.rm = TRUE)
  max.val <- max(val.mat, na.rm = TRUE)

  if(plot.individual){
    extra.bins <- c(0, 50, 100)
    plot.new()
    plot.window(xlim = c(min.x, max.x), ylim = c(1, nrow(val.mat)))
    for(i in 1:nrow(val.mat)){
        u_vals <- unique(val.mat[i,])
        u_vals <- u_vals[which(is.finite(u_vals))]
        if(length(u_vals) >= 1 && !is.na(u_vals)){
          if(length(u_vals) > 2){ 
            val.cols <- colors.from.values(val.mat[i,], use.pheatmap.colors = TRUE,
            global.color.scale = TRUE, global.min = min.val, global.max = max.val)
          }else{
            the.vals <- val.mat[i,]
            all.cols <- colors.from.values(c(the.vals, extra.bins), use.pheatmap.colors = TRUE,
            global.color.scale = TRUE, global.min = min.val, global.max = max.val)
            val.cols <- all.cols[1:length(the.vals)]
          }
          points(x = as.numeric(names(val.mat[i,])), y = rep(i, length(val.mat[i,])),
          col = val.cols, pch = 16, cex = 0.5)
        }
    }
    abline(v = 0)
    axis(1)
  }

  sum.fun <- match.fun(plot.fun)
  methyl.means <- apply(val.mat, 2, function(x) sum.fun(x, na.rm = TRUE))

  if(show.var == "se"){
    methyl.var <- apply(val.mat, 2, function(x) sd(x, na.rm = TRUE)/sqrt(length(which(!is.na(x)))))
  }

  if(show.var == "sd"){
    methyl.var <- apply(val.mat, 2, function(x) sd(x, na.rm = TRUE))
  }

  if(show.var == "var"){
    methyl.var <- apply(val.mat, 2, function(x) var(x, na.rm = TRUE))
  }

  if(plot.line){
    plot.new()
    plot.window(xlim = c(min.x, max.x), ylim = ylim)
    pos <- as.numeric(colnames(val.mat))
    plot.poly.xy(pos, methyl.means+methyl.var, pos, methyl.means-methyl.var,
    col = "#9ecae1", border = NA, new.plot = FALSE)
    points(as.numeric(colnames(val.mat)), methyl.means, type = "l")
    axis(1);axis(2)
    mtext(ylab, side = 2, line = 2)
    mtext(paste0(plot.label, "\n", nrow(val.mat), " Genes"), side = 3, line = 2)
    mtext("Relative Genomic Position", side = 1, line = 2)
  }

  if(plot.hex){
    plot.new()
    plot.window(xlim = c(min.x, max.x), ylim = ylim)
    x <- rep(as.numeric(colnames(val.mat)), nrow(val.mat))
    y <- as.vector(val.mat)
    test <- try(plot.hexbin(x,y, xlab = "Relative Genomic Position", ylab = "Percent Methylation",
    main = paste0(plot.label, "\n", nrow(val.mat),  " Genes")), silent = TRUE)
    if(class(test) == "try-error"){
      plot.text("Not enough examples to plot.")
    }
  }

}

