#This function plots values that have been centered around
#a feature by center.on.feature() The input to this function
#is a list of vectors that have all been centered by 
#center.on.feature().
#min.representation is the minimum number of genes with 
#methylation at that position that will be plotted
#seq.by indicates how finely to bin relative positions.
#leave as 1 if using bp, but change to 0.1 or 0.01 for
#relative bp
#merge.by is a smoothing parameter. After the bulk.mat
#is made, merge by combines adjacent bins so that there
#are fewer gaps in the data.
#if return.means is FALSE, the full matrix will be
#returned.

plot.centered.vals <- function(val.list, min.representation = 10, plot.label = "Center",
plot.individual = FALSE, ylim = c(0, 50), ylab = "Methylation Percent Mean",
seq.by = 1, merge.by = 1, plot.line = TRUE, plot.hex = FALSE, 
min.upstream = -2, max.downstream = 2, return.means = TRUE, verbose = FALSE){
  
  not.null <- which(sapply(val.list, function(x) !all(is.na(x))))
  if(length(not.null) == 0){stop("No values")}
  sub.val.list <- val.list[not.null]

  all.pos <- unlist(sapply(sub.val.list, function(x) as.numeric(names(x))))
  min.x <- max(c(floor(min(all.pos)), min.upstream))
  max.x <- min(c(ceiling(max(all.pos)), max.downstream))
  all.vals <- unlist(sub.val.list)
  min.val <- min(all.vals, na.rm = TRUE)
  max.val <- max(all.vals, na.rm = TRUE)

  if(plot.individual){
    extra.bins <- c(0, 50, 100)
    plot.new()
    plot.window(xlim = c(min.x, max.x), ylim = c(1, length(sub.val.list)))
    for(i in 1:length(sub.val.list)){
        u_vals <- unique(sub.val.list[[i]])
        u_vals <- u_vals[which(is.finite(u_vals))]
        if(length(u_vals) >= 1 && !is.na(u_vals)){
          if(length(u_vals) > 2){ 
            val.cols <- colors.from.values(sub.val.list[[i]], use.pheatmap.colors = TRUE,
            global.color.scale = TRUE, global.min = min.val, global.max = max.val)
          }else{
            the.vals <- sub.val.list[[i]]
            all.cols <- colors.from.values(c(the.vals, extra.bins), use.pheatmap.colors = TRUE,
            global.color.scale = TRUE, global.min = min.val, global.max = max.val)
            val.cols <- all.cols[1:length(the.vals)]
          }
          points(x = as.numeric(names(sub.val.list[[i]])), y = rep(i, length(sub.val.list[[i]])),
          col = val.cols, pch = 16, cex = 0.5)
        }
    }
    abline(v = 0)
    axis(1)
  }

  #also make a big matrix, so we can look at bulk methylation properties
  sig.fig <- length(strsplit(as.character(seq.by), "")[[1]])
  if(sig.fig > 1){sig.fig = sig.fig - 2}
  pos.seq <- seq(min.x, max.x, seq.by)
  bulk.mat <- matrix(NA, ncol = length(pos.seq), nrow = length(sub.val.list))
  rownames(bulk.mat) <- names(sub.val.list)
  colnames(bulk.mat) <- pos.seq
  for(i in 1:length(sub.val.list)){
    if(verbose){report.progress(i, length(sub.val.list))}
    rounded.pos <- round(as.numeric(names(sub.val.list[[i]])), sig.fig)
    within.bounds <- intersect(which(rounded.pos >= min.x), which(rounded.pos <= max.x))
    if(length(within.bounds) > 0){
      val.idx <- sapply(rounded.pos[within.bounds], function(x) get.nearest.pt(pos.seq, x))
      bulk.mat[i,val.idx] <- sub.val.list[[i]][within.bounds]
    }
  }

  keep.rows <- which(apply(bulk.mat, 1, function(x) !all(is.na(x))))
  bulk.mat <- bulk.mat[keep.rows,]
  #pheatmap(bulk.mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = FALSE)

  #smooth across columns to put more examples in each 
  if(merge.by > 1){
    col.list <- bin.sequence(1:ncol(bulk.mat), round(ncol(bulk.mat)/merge.by))
    merged.mat <- sapply(col.list, function(x) rowMeans(bulk.mat[,x], na.rm = TRUE))
    merged.pos <- sapply(col.list, function(x) mean(as.numeric(colnames(bulk.mat)[x])))
    colnames(merged.mat) <- merged.pos
  }else{
    merged.mat <- bulk.mat
  }

  num.examples <- apply(merged.mat, 2, function(x) length(which(!is.na(x))))
  high.rep <- which(num.examples >= min.representation)
  if(length(high.rep) == 0){stop("No genes with this high representation.")}
  sub.bulk.mat <- merged.mat[,high.rep,drop=FALSE]

  methyl.means <- colMeans(sub.bulk.mat, na.rm = TRUE)
  methyl.var <- apply(sub.bulk.mat, 2, function(x) sd(x, na.rm = TRUE)/sqrt(length(which(!is.na(x)))))

  if(plot.line){
    plot.new()
    plot.window(xlim = c(min.x, max.x), ylim = ylim)
    pos <- as.numeric(colnames(sub.bulk.mat))
    plot.poly.xy(pos, methyl.means+methyl.var, pos, methyl.means-methyl.var,
    col = "#9ecae1", border = NA, new.plot = FALSE)
    points(as.numeric(colnames(sub.bulk.mat)), methyl.means, type = "l")
    axis(1);axis(2)
    mtext(ylab, side = 2, line = 2)
    mtext(paste0(plot.label, "\n", nrow(sub.bulk.mat), " Genes"), side = 3, line = 2)
    mtext("Relative Genomic Position", side = 1, line = 2)
  }

  if(plot.hex){
    plot.new()
    plot.window(xlim = c(min.x, max.x), ylim = ylim)
    x <- rep(as.numeric(colnames(sub.bulk.mat)), nrow(sub.bulk.mat))
    y <- as.vector(sub.bulk.mat)
    test <- try(plot.hexbin(x,y, xlab = "Relative Genomic Position", ylab = "Percent Methylation",
    main = paste0(plot.label, "\n", nrow(sub.bulk.mat),  " Genes")), silent = TRUE)
    if(class(test) == "try-error"){
      plot.text("Not enough examples to plot.")
    }
  }

  if(return.means){
    names(methyl.means) <- colnames(sub.bulk.mat)
    invisible(methyl.means)
  }else{
    invisible(sub.bulk.mat)
  }

}

