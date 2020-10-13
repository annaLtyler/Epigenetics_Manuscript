#This function compares different methylation matrices from around
#a single gene, for example in treatment and control animals.
#it compares only the common positions between the two matrices

compare.methyl.results <- function(methyl.results1, methyl.results2, 
methyl.labels = c("Methyl1", "Methyl2"), min.change = 50, 
save.images = FALSE, path = ".", plot.label = ""){

    mat1 <- methyl.results1$methylation.matrix
    mat2 <- methyl.results2$methylation.matrix

    common.pos <- intersect(colnames(mat1), colnames(mat2))

    mat1.idx <- match(common.pos, colnames(mat1))
    mat2.idx <- match(common.pos, colnames(mat2))

    common.ind <- intersect(rownames(mat1), rownames(mat2))
    mat1.ind.locale <- match(common.ind, rownames(mat1))
    mat2.ind.locale <- match(common.ind, rownames(mat2))
    
    mat1 <- mat1[mat1.ind.locale,mat1.idx]
    mat2 <- mat2[mat2.ind.locale,mat2.idx]

    diff.mat <- mat1 - mat2
    diff.mat[which(diff.mat == 0)] <- NA

    if(save.images){
        jpeg(file.path(path, paste0("Methylation.Differences.", plot.label, ".jpg")),
        width = 10, height = 5, units = "in", res = 300)
    }else{
        quartz(width = 10, height = 5)
    }
 
    strain.cols <- t(apply(diff.mat, 1, 
    function(x) colors.from.values(x, use.pheatmap.colors = TRUE, global.color.scale = TRUE, global.min = -100, global.max = 100)))
    methyl.pos <- as.numeric(colnames(diff.mat))

    layout(matrix(c(1,2), ncol = 1), heights = c(1, 0.3))
    plot.new()
    plot.window(xlim = c(min(methyl.pos), max(methyl.pos)), 
    ylim = c(0,nrow(diff.mat)+1))
    par(mar = c(0,4,2,4))
    plot.dim <- par("usr")
    plot.width = plot.dim[2] - plot.dim[1]
    par(xpd = TRUE)
    for(i in 1:nrow(diff.mat)){
        points(x = methyl.pos, y = rep(i, length(methyl.pos)), 
        pch = "|", col = strain.cols[i,], cex = 1.5)
        text(plot.dim[2], i, rownames(diff.mat)[i])
    }
    par(xpd = FALSE)
    mtext("Methylation", side = 2)

    par(mar = c(2,2,2,2))
    min.meth.change = min(diff.mat, na.rm = TRUE)
    max.meth.change = max(diff.mat, na.rm = TRUE)
    imageWithTextColorbar(matrix(c(min.meth.change:max.meth.change), ncol = 1), use.pheatmap.colors = TRUE, 
    cex = 1, orientation = "h")

    if(save.images){dev.off()}
   
    expr1 <- methyl.results1$expr
    expr2 <- methyl.results2$expr
        
    expr.ind.locale1 <- match(common.ind, colnames(expr1))
    expr.ind.locale2 <- match(common.ind, colnames(expr2))

    expr1 <- expr1[,expr.ind.locale1]
    expr2 <- expr2[,expr.ind.locale2]
    
    if(save.images){
        jpeg(file.path(path, paste0("Expression.Differences.", plot.label, ".jpg")),
        width = 6, height = 4, units = "in", res = 300)
    }else{
        quartz(width = 6, height = 4)
    }
    plot.grouped.boxes(list(expr1, expr2), type = "matrix", group.labels = methyl.labels)
    if(save.images){dev.off()}

    #is there a relationship between change in methylation 
    #and the change in expression?
    test.locale <- which(abs(diff.mat) >= min.change, arr.ind = TRUE)
    beta.tests <- matrix(NA, nrow = nrow(test.locale), ncol = 3)
    colnames(beta.tests) <- c("strain", "position", "beta")
    for(i in 1:nrow(beta.tests)){
        strain.idx <- test.locale[i,1]
        position.idx <- test.locale[i,2] 
        strain.expr1 <- expr1[,strain.idx]
        strain.expr2 <- expr2[,strain.idx]
        beta.tests[i,1] <- colnames(expr1)[strain.idx]
        #boxplot(list(strain.expr1, strain.expr2))
        exprv <- c(strain.expr1, strain.expr2)
        
        methyl1 <- mat1[strain.idx,position.idx]
        methyl2 <- mat2[strain.idx,position.idx]
        beta.tests[i,2] <- colnames(mat1)[position.idx]
        
        dummy.methyl1 <- rep(methyl1, length(strain.expr1))
        dummy.methyl2 <- rep(methyl2, length(strain.expr2))
        methylv <- c(dummy.methyl1, dummy.methyl2)
        #plot(exprv, methylv, pch = 16)
        model <- lm(exprv~methylv)
        beta.tests[i,3] <- as.numeric(coef(model)[2])
    }
    u_strain <- unique(beta.tests[,1])
    beta.by.strain <- lapply(u_strain, function(x) beta.tests[which(beta.tests[,1] == x),2:3])    
    names(beta.by.strain) <- u_strain    

    if(save.images){
        jpeg(file.path(path, paste0("Methylation.Change.Effects.", plot.label, ".jpg")),
        width = 5, height = 5, units = "in", res = 300)
        }else{
            quartz()
        }
    beta.val <- lapply(beta.by.strain, function(x) as.numeric(x[,2]))

    stripchart(beta.val, pch = 16, vertical = TRUE, method = "stack", col = "gray",
    offset = 0.5)
    abline(h = 0)
    if(save.images){dev.off()}
}

