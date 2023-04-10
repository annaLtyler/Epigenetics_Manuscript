#This function tests different dimension-reduction
#methods on a given data matrix
#The data matrix should have observations 
#(e.g. genes) in rows and independent 
#samples (e.g. mice) in columns.

test_dim_red <- function(data.mat, col = "black", pdf.filename = "~/Desktop/test.pdf"){

    require(cluster);require(umap);require(fastICA);require(tsne)
    require(vegan);require(Rdimtools)

    pdf(pdf.filename, width = 8, height = 8)
    par(mfrow = c(2,2))
    #trying out different low-dimenstional visualizations
    data.pc <- plot.decomp(t(data.mat), plot.results = FALSE, pc = 4)
    data.decomp <- data.pc$u
    plot(data.decomp[,c(1,2)], pch = 16, col = col, xlab = "PC1", ylab = "PC2",
    main = "SVD of data matrix")

    data.cor <- cor(data.mat)
    plot.decomp(data.cor, cols = col, main = "SVD of correlation matrix")

    umap.defaults$n_neighbors = 10
    umap.defaults$min_dist = 0.01
    umap.defaults$spread = 2
    data.umap <- umap(data.cor)
    data.decomp <- data.umap$layout
    plot(data.decomp, pch = 16, col = col, xlab = "UMAP1", ylab = "UMAP2",
        main = "UMAP")

    data.ica <- fastICA(data.mat, n.comp = 2)
    plot(data.ica$K, col = col, pch = 16, main = "ICA")

    data.tsne <- tsne(dist(t(data.mat)), k = 2, initial_dims = 9, perplexity = 2)
    plot(data.tsne, col = col, pch = 16, xlab = "TSNE1", ylab = "TSNE2",
        main = "tSNE")


    data.mds <- cmdscale(dist(t(data.mat)), 2)
    plot(data.mds, col = col, pch = 16, main = "MDS")


    data.iso <- isomap(dist(t(data.mat)), k = 9)
    plot(data.iso$points[,1:2], pch = 16, col = col, main = "isomap")

    
    data.rnd <- do.rndproj(t(data.mat), ndim = 2)
    plot(data.rnd$Y, pch = 16, col = col, main = "Random Projection")

    data.lle <- do.lle(t(data.mat), ndim = 2)
    plot(data.lle$Y, pch = 16, col = col, main = "Local Linear Embedding")

    dev.off()
}