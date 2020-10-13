#This function returns a single methylated island close to either
#the TSS or the TES. For the TSS, the algorithm does the following steps:
#1. find island overlapping TSS
#2. if none, look at the average methylation in all islands.
#3. If any have an average methylation less than tss.min.methyl, it
# takes the one closest to the TSS. 
# If there are multiple close islands, the largest one is taken.
# If there are no hypomethylated islands, it takes the island
# nearest the TSS regardless of the average methylation
# Types of islands selected are "overlap",  "nearest-large-low-island",
# "nearest-low-island" and "nearest-island"
#If we are looking for the island near the TES, and there isn't
#one overlapping, we just take the next one upstream regardless of 
#average methylation.
# types of islands selected are "overlap" or "upstream"

pick.island <- function(gene.name, island.obj, gene.info, anchor = c("TSS", "TES"), 
tss.min.methyl = 10){

    if(length(island.obj) < 2){
        return(NA)
    }

    gene.locale <- which(gene.info[,"external_gene_name"] == gene.name)
    gene.start <- unique(gene.info[gene.locale,"start_position"])[1]
    gene.end <- unique(gene.info[gene.locale,"end_position"])[1]
    gene.strand <- unique(gene.info[gene.locale,"strand"])
    pos <- as.numeric(names(island.obj[[1]]))
    
    if(gene.strand == 1){
        if(anchor == "TSS"){
            nearest.pos <- get.nearest.pt(gene.start, pos)
        }
        if(anchor == "TES"){
            nearest.pos <- get.nearest.pt(gene.end, pos)
        }
    }

    if(gene.strand == -1){
        if(anchor == "TSS"){
            nearest.pos <- get.nearest.pt(gene.end, pos)
        }
        if(anchor == "TES"){
            nearest.pos <- get.nearest.pt(gene.start, pos)
        }
    }

    island.id <- island.obj[[1]][nearest.pos] 
    all.island.idx <- which(island.obj[[1]] == island.id)

    island.pos <- pos[all.island.idx]
    if(gene.strand == 1){
        if(anchor == "TSS"){
            does.overlap <- (gene.start > min(island.pos) && gene.start < max(island.pos))
        }
        if(anchor == "TES"){
            does.overlap <- (gene.end > min(island.pos) && gene.end < max(island.pos))
        }
    }

    if(gene.strand == -1){
        if(anchor == "TSS"){
            does.overlap <- (gene.end > min(island.pos) && gene.end < max(island.pos))
            #plot(island.pos);abline(h = gene.end)
        }
        if(anchor == "TES"){
            does.overlap <- (gene.start > min(island.pos) && gene.start < max(island.pos))
            #plot(island.pos);abline(h = gene.end)
        }
    }

    #if the island does overlap the region of interest, pick this one.
    if(does.overlap){
        selected.island <- list(island.obj[[1]][all.island.idx], island.obj[[2]][[island.id]], "overlap")
    }else{
        #if the island does not overlap the region of interest
        #look for hypomethylated islands if the anchor is the TSS
        if(anchor == "TSS"){
            #find the average methylation values of all islands
            island.means <- sapply(island.obj[[2]], function(x) mean(x, na.rm = TRUE))
            #barplot(island.means)
            island.vals <- island.obj[[1]]
            u_islands <- unique(island.vals)
            island.size <- sapply(u_islands, function(x) length(which(island.vals == x)))
            island.min.pos <- sapply(u_islands, function(x) min(pos[which(island.vals == x)]))
            island.max.pos <- sapply(u_islands, function(x) max(pos[which(island.vals == x)]))
            low.islands <- which(island.means <= tss.min.methyl)
            if(length(low.islands) > 0){ #if there are some hypomethylated islands
    
                if(gene.strand == 1){
                    gene.anchor <- gene.start
                }else{
                    gene.anchor <- gene.end
                }
                nearest.island.min <- get.nearest.pt(island.min.pos[low.islands], gene.anchor)
                nearest.island.max <- get.nearest.pt(island.max.pos[low.islands], gene.anchor)
                nearest.island <- unique(c(nearest.island.min, nearest.island.max))
                if(length(nearest.island) > 1){
                    #pick the largest of the two nearest islands.
                    large.island.locale <- which.max(island.size[nearest.island])
                    island.id <- low.islands[nearest.island[large.island.locale]]
                    tag <- "nearest-large-low-island"
                }else{
                    island.id <- low.islands[nearest.island]
                    tag <- "nearest-low-island"
                }
            }else{
                #if there are no hypomethylated islands, just keep the first one we picked
                island.id <- island.id
                tag <- "nearest-island"
            }
        }else{
            #If we are looking for the TES instead, just get the island
            #next upstream of the TES.
            if(gene.strand == 1){
                island.id <- island.id - 1
            }
            if(gene.strand == -1){
                island.id <- island.id + 1
            }
            tag <- "upstream"
        }

    if(island.id == 0 || island.id > max(island.obj[[1]])){
        return(NA)
    }

        all.island.idx <- which(island.obj[[1]] == island.id)    
        selected.island <- list(island.obj[[1]][all.island.idx], island.obj[[2]][[island.id]], tag)
    }

    return(selected.island)

}