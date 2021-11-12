#This function takes two bed-formated tables and 
#uses the enrichment function from ChromHMM to 
#determine whether the elements of the first bed-formated 
#tables are enriched based on the positions in the second
#bed-formated table. The formula is as stated in the ChromHMM
#manual: https://usermanual.wiki/Document/ChromHMMmanual.899564556/help
#By default the fold enrichment calculation is as follows, let:
#A - the number of bases in the state
#B - the number of bases in the external annotation
#C - the number of bases in the state and the external annotation
#D - the number of bases in the genome
#The fold enrichment is then defined as (C/A)/(B/D).
#test.bed <- "~/Documents/Projects/Epigenetics/Epigenetics_Manuscript/Data/mm10/CTCF.bed.gz"
#bed.file2 <- "~/Documents/Projects/Epigenetics/Epigenetics_Manuscript/Data/mm10/enhancers.bed.gz"
#chrom.size.file <- "/Users/atyler/Documents/Projects/Epigenetics/Epigenetics_Manuscript/Data/support_files/mm10.txt"
enrichment_bed <- function(state.bed, annotation.bed, chrom.size.file){
    require(data.table)

    bedz1 = gzfile(state.bed,'rt')  
    bed1 = read.table(bedz1, header = F, stringsAsFactors = FALSE, fill = TRUE)

    bedz2 = gzfile(annotation.bed,'rt')  
    bed2 = read.table(bedz2, header = F, stringsAsFactors = FALSE, fill = TRUE)

    closeAllConnections()

    chr <- unique(bed1[,1])

    chrom.size <- read.table(chrom.size.file, header = F, stringsAsFactors = FALSE)
    chr.locale <- match(chr, chrom.size[,1])

    A <- sum(bed1[,3] - bed1[,2]+1) #A - the number of bases in the state
    B <- sum(bed2[,3] - bed2[,2]+1) #B - the number of bases in the external annotation
    D <- sum(chrom.size[chr.locale,2]) #D - the number of bases in the genome
    
    #To calculate C, we need to calculate all pair-wise overlaps
    #of the coordinates and count up all the overlapping bases.
    dt1 = data.table(chrom = bed1[,1], start = bed1[,2], end = bed1[,3])
    dt2 = data.table(chrom = bed2[,1], start = bed2[,2], end = bed2[,3])
    setkey(dt1, chrom, start, end)
    setkey(dt2, chrom, start, end)
    
    overlap.table <- foverlaps(dt1, dt2, type="any")
    overlap.only <- overlap.table[which(!is.na(overlap.table[,3])),]

    get_num_overlap <- function(overlap.table.row, plot.results = FALSE){
        region1 <- as.numeric(overlap.table.row[2:3])
        region2 <- as.numeric(overlap.table.row[4:5])
        sorted.coord <- sort(c(region1, region2))
        overlap.size <- sorted.coord[3] - sorted.coord[2]
        
        if(plot.results){
            quartz()
            plot.new()
            plot.window(xlim = c(min(sorted.coord), max(sorted.coord)), ylim = c(0,2))
            segments(x0 = region1[1], x1 = region1[2], y0 = 1)
            segments(x0 = region2[1], x1 = region2[2], y0 = 1.5)
            abline(v = c(sorted.coord[2:3]))
        }
        return(overlap.size)
    }

    all.overlap <- apply(overlap.only, 1, get_num_overlap)
    C <- sum(all.overlap)

    enrichment.score <- (C/A)/(B/D)
    return(enrichment.score)
}
