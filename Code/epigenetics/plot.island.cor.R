#This function takes in results from methyl.expression.cor
#and plots the relationship between a single specified
#island and the expression of a gene

plot.island.cor <- function(methyl.results, island.which = 1, col.table){

    island.locale <- which(methyl.results$island.bins == island.which)
    sub.methyl <- methyl.results$methylation.matrix[,island.locale,drop=FALSE]
    methyl.means <- rowMeans(sub.methyl, na.rm = TRUE)       

    expr.order <- order.strains(names(methyl.means), colnames(methyl.results$expr), col.table)
    ordered.expr <- methyl.results$expr[,expr.order]
    cols <- col.table[order.strains(names(methyl.means), col.table[,1], col.table),3]

    plot.new()
    plot.window(xlim = c(min(methyl.means, na.rm = TRUE), max(methyl.means, na.rm = TRUE)),
    ylim = c(min(ordered.expr), max(ordered.expr)))
    for(i in 1:length(methyl.means)){
        points(x = rep(methyl.means[i], nrow(ordered.expr)), y = ordered.expr[,i],
        col = cols[i], pch = 16, cex = 1.5)
    }
    axis(1);axis(2)
    mtext("Percent Methylated", side = 1, line = 2)
    mtext("Expression", side = 2, line = 2)

}