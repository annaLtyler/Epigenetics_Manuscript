#This function takes in EdgeR output (qlf object)
#It filters the results by p value and log fold 
#change, and adds additional about the genes using 
#biomaRt. The results are written to a file for 
#further processing

EdgeR.to.table <- function(qlf, result.name, lFC = 0.5, pval = 0.05){
    
  #Convert list into table
  subtable <- qlf$table
  if(!is.null(lFC)){
  	subtable <- subtable[which(abs(subtable$logFC) >= lFC),]
  	}
  if(!is.null(pval)){
	subtable <- subtable[which(subtable$PValue <= pval),]
  	}
  
  #Get the ensembl ids for the de genes
  ensIds=rownames(subtable)
  
  # get chromosomal coordinates for the mapped genes 
  attr=c('ensembl_gene_id','external_gene_name', 'chromosome_name','start_position','end_position')
  chrlocs=getBM(attributes=attr, filters = "ensembl_gene_id", values=ensIds, mart=mm10)
  rownames(chrlocs)=chrlocs$ensembl_gene_id
  chrlocs=chrlocs[,-1]
  
  # merge this with the de results table
  RNA = merge.data.frame(chrlocs, subtable, by.x=0, by.y=0, all.x=TRUE, all.y=TRUE)
  colnames(RNA) = c("EnsemblID","Gene", "Chr","Start","End", "logFC", "logCPM", "F", "PValue")
  # Remove NAs
  RNA = as.data.frame(RNA[!is.na(RNA$Chr), ])
  #Remove extra column w/ row numbers
  RNA <- RNA["Chr"!= "Y"]
  RNA <- RNA["Chr"!= "MT"]
  
  write.table(RNA, file = result.name, row.names = FALSE, sep = "\t", quote = FALSE)

}