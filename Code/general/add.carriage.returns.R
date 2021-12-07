#This function splits up a text message by words 
#making sure that there are no more than the specified
#number of characters on any line.
#This function assumes words are separated by spaces.

add.carriage.returns <- function(the.text, max.char.per.line = 40){

    words <- strsplit(the.text, " ")[[1]]
    word.len <- sapply(words, nchar)

    if(sum(word.len) <= max.char.per.line){
        return(the.text)
    }

    split.line <- list()
    running.total <- cumsum(word.len)
    start.idx <- 1
    sub.words <- words
    while(length(running.total) > 0){
        this.line <- which(running.total <= max.char.per.line)
        split.line[[start.idx]] <- paste0(paste(sub.words[this.line], collapse = " "), "\n")
        the.rest <- setdiff(1:length(running.total), this.line)
        sub.words <- sub.words[the.rest]
        word.len <- sapply(sub.words, nchar)
        running.total <- cumsum(word.len)
        start.idx <- start.idx + 1
    }

    final.line <- Reduce(paste, split.line)
    return(final.line)

}