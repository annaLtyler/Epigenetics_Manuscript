get_one_geno  <- function(gene.name, transcript.info, transcript.haplotypes, 
chrom.states, strain.key, geno_type = c("chromatin", "genotype"), 
perm.order = NULL){
  
  gene.id <- transcript.info[which(transcript.info[,"external_gene_name"] == gene.name),"ensembl_gene_id"][1]

  geno_type <- geno_type[1]
  info.locale <- which(transcript.info[,"ensembl_gene_id"] == gene.id)[1]

  if(is.na(info.locale)){return(NA)}

  gene.chr <- transcript.info[info.locale,"chromosome_name"]
  gene.start <- transcript.info[info.locale,"start_position"]

  if(geno_type == "chromatin"){
    one.geno <- list(calc.chrom(gene.name, transcript.info, transcript.haplotypes, 
    chrom.states, strain.key = strain.key, perm.order = perm.order))
  }else{
    gene.locale <- which(names(transcript.haplotypes) == gene.id)
    if(length(gene.locale) == 1){return(NA)}
    one.geno <- list(transcript.haplotypes[[gene.locale]])
  }

  #add attributes
    names(one.geno) <- gene.chr
    if(gene.chr == "X"){
      attr(one.geno, "is_x_chr") <- TRUE
    }else{
      attr(one.geno, "is_x_chr") <- FALSE
    }
    attr(one.geno, "crosstype") <- "do"
    attr(one.geno, "alleles") <- LETTERS[1:8]
    attr(one.geno, "alleleprobs") <- TRUE
    attr(one.geno, "class") <- c("calc_genoprob", "list")
    return(one.geno)
}
