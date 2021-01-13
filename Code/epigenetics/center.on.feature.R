#center a result on a given feature and orient
#genes such that all transcription is going
#in the same direction
#this function expects that the names of vals are the
#genomic coordinates
#if feature == "full", the upstream, downstream, and
#gene body coordinates are set as negative, above one,
#and from 0 to one respectively.

center.on.feature <- function(gene.name, gene.info, vals, 
feature = c("tss", "tes", "1stexon", "promoter", "full")){

  if(is.na(gene.name)){return(NA)}

  regions <- c("1stexon", "promoter", "full")
  is.region <- as.logical(length(which(regions == feature)))

  if(length(vals) == 1){
    return(vals)
  }

  if(is.null(names(vals))){names(vals) <- vals}

  gene.locale <- which(gene.info[,"external_gene_name"] == gene.name)
  if(length(gene.locale) == 0){return(NA)}
  gene_strand <- unique(gene.info[gene.locale,"strand"])
  pos <- as.numeric(names(vals))

  #plot(pos, vals)
  #tss <- unique(gene.info[gene.locale,"start_position"])
  #tes <- unique(gene.info[gene.locale,"end_position"])
  #abline(v = c(tss, tes))

  if(feature == "tss"){
    feature.end = NA
    if(gene_strand == 1){
      feature.start <- unique(gene.info[gene.locale,"start_position"])
    }else{
      feature.start <- unique(gene.info[gene.locale,"end_position"])
    }
  }

  if(feature == "tes"){
    feature.end <- NA
    if(gene_strand == 1){
      feature.start <- unique(gene.info[gene.locale,"end_position"])
    }else{
      feature.start <- unique(gene.info[gene.locale,"start_position"])
    }
  }

  if(feature == "1stexon"){
    if(gene_strand == 1){
      #first.locale = 1
      first.locale <- which.min(gene.info[gene.locale,"exon_chrom_start"])
      feature.start <- gene.info[gene.locale[first.locale],"exon_chrom_start"]
      feature.end <- gene.info[gene.locale[first.locale],"exon_chrom_end"]
    }else{
      #first.locale <- nrow(gene.info[gene.locale,])
      first.locale <- which.max(gene.info[gene.locale,"exon_chrom_end"])
      feature.start <- gene.info[gene.locale[first.locale],"exon_chrom_end"]
      feature.end <- gene.info[gene.locale[first.locale],"exon_chrom_start"]
    }
  }

if(feature == "promoter"){
    if(gene_strand == 1){
      tss <- gene.info[gene.locale[1],"start_position"]
      feature.start <- tss - 700
      feature.end <- tss + 200
    }else{
      tss <- gene.info[gene.locale[1],"end_position"]
      feature.start <- tss - 200
      feature.end <- tss + 700
    }
  }

   if(feature == "full"){
    if(gene_strand == 1){
      feature.start <- gene.info[gene.locale[1],"start_position"]
      feature.end <- gene.info[gene.locale[1],"end_position"]
    }else{
      feature.start <- gene.info[gene.locale[1],"end_position"]
      feature.end <- gene.info[gene.locale[1],"start_position"]
    }
  }

  #abline(v = c(feature.start, feature.end), col = c("red", "blue"))

    #flip if we are on the negative strand
    if(gene_strand == 1){
      rel.pos <- pos - feature.start
    }else{
      rel.pos <- -pos + feature.start
    }
  

  if(is.region){
    #scale by the length of the region
    scaled.pos <- rel.pos/(abs(feature.end - feature.start)) 
    names(vals) <- scaled.pos
  }else{
    names(vals) <- rel.pos
  }

  return(vals)
}