#This function reads in a database file format from 
#the Broad for doing GSEA. These files list genes 
#associated with GO terms, but sometimes the number of
#genes is too great for the row, so they spill over into
#the next line. This function parses those spillovers 
#to associate genes with each GO term. It returns a
#list. Each element contains the genes for one GO term.

read.gmt <- function(filename){

    pre.table <- as.matrix(read.table(filename, sep = "\t", stringsAsFactors = FALSE, fill = TRUE))
    
    term.col <- pre.table[,2]
    term.locale <- grep("GO:", term.col)
    term.ids <- pre.table[term.locale,1]

    all.go.genes <- vector(mode = "list", length = length(term.locale))
    names(all.go.genes) <- term.ids

    for(i in 1:length(term.locale)){
        go.start <- term.locale[i]
        if(i == length(term.locale)){
            go.end <- go.start
        }else{
            go.end <- term.locale[(i+1)] - 1
        }
        if(go.start == go.end){
            go.genes <- pre.table[go.start,3:ncol(pre.table)]
            go.genes <- go.genes[which(go.genes != "")]
        }else{
            go.genes <- pre.table[go.start,3:ncol(pre.table)]
            go.start <- go.start + 1
            add.genes <- as.vector(pre.table[go.start:go.end,])
            go.genes <- c(go.genes, add.genes)
            go.genes <- go.genes[which(go.genes != "")]
        }
        all.go.genes[[i]] <- go.genes
    }

    return(all.go.genes)

}