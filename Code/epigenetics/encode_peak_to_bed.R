#This function reads in a set of peak position files from mouse encode 
#http://chromosome.sdsc.edu/mouse/download.html
#and converts the positions to bed files that we 
#can use to look for positional enrichments of other
#genomic features.
#source_dir <- "~/Downloads/liver"
#target_dir <- here("Data", "mm10")
encode_peak_to_bed <- function(source_dir, target_dir, bp.buffer = 5){
    
    peak.files <- list.files(source_dir, full.names = TRUE)
    for(i in 1:length(peak.files)){
        peak.table <- read.table(peak.files[i], stringsAsFactors = FALSE, header = FALSE)
        #add buffer
        start.pos <- peak.table[,2] - bp.buffer
        end.pos <- peak.table[,2] + bp.buffer
        chr <- unique(peak.table[,1])
        bed.table <- cbind(peak.table[,1], start.pos, end.pos)
        bed.filename <- gsub("txt", "bed", basename(peak.files[i]))
        write.table(bed.table, file.path(target_dir, bed.filename), quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
        system(paste("gzip", file.path(target_dir, bed.filename)))
        }
}