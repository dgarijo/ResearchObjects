library(customProDB)


calcRPKMValues <- function(){
  args <- commandArgs(TRUE)
  bam.file <- args[1]
  cpdb.file <- args[2]
  fastq.file <- args[3]
  annotation.path <- args[4]
  output.file <- args[5]
  prot.colname <- args[6]
  sampleid.colname <- args[7]
  rpkm.colname <- args[8]
  
  load(paste(annotation.path, "/exon_anno.RData", sep = ""))
  load(paste(annotation.path, "/ids.RData", sep = ""))
  
  if(is.na(prot.colname)) {
    prot.colname <- "proname"
  }
  
  if(is.na(sampleid.colname)) {
    sampleid.colname <- "sample_id"
  }
  
  if (is.na(rpkm.colname)) {
    rpkm.colname <- "RPKM"
  }
  
  #Brute force -> correct in future
  sample.id <- gsub("(.*)\\..*", "\\1", basename(fastq.file))
  cpdb.df <- read.delim(cpdb.file, header = T, stringsAsFactors = FALSE)
  
  # calculate RPKM values from a BAM file based on exon read counts.
  rpkm.values <- calculateRPKM(bamFile = bam.file,
                                exon = exon,
                                ids = ids)
  
  # calculateRPKM returns a vector with corresponding protein_ids as names
  exonrc.rpkm.df <- setNames(data.frame(names(rpkm.values), round(rpkm.values,2)),
                             c(prot.colname, rpkm.colname))
  cpdb.rpkm.df <- merge(cpdb.df, exonrc.rpkm.df, all.x=TRUE)
  cpdb.rpkm.df[,sampleid.colname] <- sample.id
  prot.rpkm.df <- unique(cpdb.rpkm.df[,c(prot.colname, sampleid.colname, rpkm.colname)])
  write.table(prot.rpkm.df, output.file, sep = "\t", row.names = F, quote = F)
  return(output.file)
}

calcRPKMValues()