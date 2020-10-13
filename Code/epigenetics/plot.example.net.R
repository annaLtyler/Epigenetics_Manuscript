plot.example.net <- function(Lij, gene.ids, mart, layout.call = layout_nicely){
	
	require(igraph)

	gene.names <- getBM(c("entrezgene", "external_gene_name"), "entrezgene", gene.ids, mart)
	names.by.locus <- gene.ids
	for(i in 1:length(names.by.locus)){
		gene.locale <- which(gene.names[,1] == gene.ids[i])
		if(length(gene.locale) > 0){
			names.by.locus[i] <- gene.names[gene.locale,2]	
			}
		}
	

	na.locale <- which(is.na(names.by.locus))
	names.by.locus[na.locale] <- names(names.by.locus)[na.locale]

	adj.mat <- Lij
	adj.mat[which(is.na(adj.mat))] <- 0
	
	colnames(adj.mat)[1:length(names.by.locus)] <- rownames(adj.mat)[1:length(names.by.locus)] <- names.by.locus
	net <- graph.adjacency(adj.mat, mode = "directed", weighted = TRUE)
	e.col <- rep("#d8b365", ecount(net))
	e.col[which(E(net)$weight < 0)] <- "#5ab4ac"
	
	E(net)$arrow.size <- 0.4
	E(net)$width <- 1+abs(E(net)$weight)/12

	the.layout <- layout_nicely(net)
	norm.layout <- norm_coords(the.layout, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
	stretched.layout <- norm.layout*1.2
	stretched.layout[,1] <- stretched.layout[,1] + 0.1 #move along x axis
	stretched.layout[,2] <- stretched.layout[,2] + 0.2 #move along y axis

	igf.interactors <- unique(c(names(which(adj.mat["Final.IGF.1",] != 0)), names(which(adj.mat[,"Final.IGF.1"] != 0))))
	to.mark <- c("Final.IGF.1", igf.interactors)

	vcol <- rep(adjustcolor("#7fc97f", alpha = 0.6), vcount(net))
	vcol[match(to.mark, V(net)$name)] <- adjustcolor("#beaed4", alpha = 0.6)
	V(net)$color <- vcol

	modes = c("all", "out", "in")
	for(m in 1:length(modes)){
		deg <- degree(net, mode=modes[m])
		V(net)$size <- deg*3

		quartz(width = 8, height= 8)
		par(xpd = TRUE)
		plot(net, edge.curved = 0.1, layout = stretched.layout, vertex.frame.color = "#ffffff", edge.color = e.col, rescale = FALSE, main = modes[m])
		par(xpd = FALSE)
		}


		layout.matrix <- readRDS("layout.matrix.RData")
		quartz(width = 8, height= 8)
		par(xpd = TRUE)
		plot(net, edge.curved = 0.1, layout = layout.matrix, vertex.frame.color = "#ffffff", edge.color = e.col)
		par(xpd = FALSE)
		

	
# # 	if(layout.call == "manual"){
		# tkp.id <- tkplot(net, layout = stretched.layout)
		# done <- readline(prompt = "Press return when ready:\n")
		# layout.matrix <- tkplot.getcoords(tkp.id)
		# tkplot.close(tkp.id)
		# }


saveRDS(layout.matrix, "layout.matrix.RData")
		plot(net, edge.curved = 0.1, layout = layout.matrix, vertex.frame.color = "#ffffff", edge.color = e.col)


	sex.interactors <- c("sex", unique(c(names(which(adj.mat["sex",] != 0)), names(which(adj.mat[,"sex"] != 0)))))

	layers <- rep(NA, vcount(net))
	names(layers) <- V(net)$name
	layers[sex.interactors] <- 3	
	layers[to.mark] <- 2
	na.left <- which(is.na(layers))
	layers[intersect(na.left, which(deg < 4))] <- 1
	na.left <- which(is.na(layers))
	layers[intersect(na.left, which(deg >= 4))] <- 4
	cbind(layers, deg)	
	
	test.layout <- layout_with_sugiyama(simplify(net), layers = layers)
	plot(net, edge.curved = 0.1, layout = test.layout$layout, vertex.frame.color = "#ffffff", edge.color = e.col)
	
 	a <- heatmap(Lij)	
 	Lij[which(Lij == 0)] <- NA
 	imageWithText(Lij[a$rowInd, a$colInd], show.text = FALSE)

}