#This function takes in a matrix that can be converted
#to a bipartite graph. The matrix should have 0's where
#there are no edges and a numerical value where there
#are edges. Row names and column names are used to name
#edges. The function creates the bipartite graph,
#and optionally plots one or both of the projections.

plot_bipartite <- function(bip_mat, plot.results = TRUE, vertex.label.cex = 1,
vertex.label.dist = 1, max.edge.lwd = 3){

    #find all connected pairs
    bip.pairs <- which(bip_mat != 0, arr.ind = TRUE)
    bip.weights <- bip_mat[which(bip_mat != 0)]
    if(!is.null(rownames(bip_mat))){
        bip.pairs[,1] <- rownames(bip_mat)[as.numeric(bip.pairs[,1])]
        bip.pairs[,2] <- colnames(bip_mat)[as.numeric(bip.pairs[,2])]
    }

    inc.mat <- incidence.matrix(bip.pairs)
    inc.net <- graph_from_incidence_matrix(inc.mat)
    E(inc.net)$weight <- bip.weights
    proj <- bipartite_projection(inc.net)

    if(plot.results){
        
        plot_proj <- function(net_proj, label.vertex = NULL){
            edge.strength <- E(net_proj)$weight
            edge.scale.factor <- max.edge.lwd/max(edge.strength)
            edge.lwd <- edge.strength*edge.scale.factor

            #quartz()
            if(length(unique(edge.strength)) > 1){
                edge.col <- colors.from.values(edge.strength, use.pheatmap.colors = TRUE)
            }else{
                edge.col <- "lightblue"
            }
            
            vertex.color = rep("lightblue", vcount(net_proj))
            if(!is.null(label.vertex)){
                vertex.locale <- which(V(net_proj)$name %in% label.vertex)
                if(length(vertex.locale) > 0){
                    vertex.color[vertex.locale] <- "red"
                }
            }

            #pdf("~/Desktop/net.pdf", width = 15, height = 15)
            plot(net_proj, layout = layout_nicely, vertex.size = 5, 
            vertex.color = vertex.color, edge.color = edge.col, 
            vertex.label.cex = vertex.label.cex, edge.width = edge.lwd,
            vertex.label.dist = vertex.label.dist)
            #dev.off()
        }
    
        plot_proj(inc.net, label.vertex = colnames(bip_mat))
        plot_proj(proj[[1]])
        plot_proj(proj[[2]])
    }
    
        
    result <- list("Full_Network" = inc.net, "Projection1" = proj[[1]], "Projection2" = proj[[2]])
    invisible(result)

}
