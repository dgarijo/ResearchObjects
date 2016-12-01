# library(reshape2)
# library(seqinr)
# source("/local/home/ravali/genericFunctionsPGomics.R")

getProteinLengthFromRefseqID <- function(protein.fastafile,
                                         refseq.protein.groups,
                                         prot.colname,
                                         one.prot.colname,
                                         prot.length.colname,
                                         group.proteins.max.or.min = "max") {
  split.protein.groups.list <- strsplit(refseq.protein.groups, ", ")
  names(split.protein.groups.list) <- refseq.protein.groups
  refseq.protein.ids.df <- setNames(melt(split.protein.groups.list),
                                    c("labkey_prot_name", prot.colname))
  
  # Get protein info (protein name, length, tx name, gene name) from the fasta file:
  fasta.refseq.prot.aaseqs.list <- read.fasta(file = protein.fastafile,
                                              seqtype = "AA")
  refseq.prots.names <- unlist((getName(fasta.refseq.prot.aaseqs.list)))
  refseq.prots.description <- unlist(getAnnot(fasta.refseq.prot.aaseqs.list))
  refseq.prots.length <- unlist(getLength(fasta.refseq.prot.aaseqs.list))
  refseq.id.prot.info.df <- setNames(data.frame(refseq.prots.names,
                                                refseq.prots.description,
                                                refseq.prots.length,
                                                stringsAsFactors = FALSE),
                                     c(one.prot.colname, prot.annot.colname,
                                       prot.length.colname))
  
  refseq.id.length.df <- refseq.id.prot.info.df
  # Partial match refseq ID # Enable-ing partial search for fasta headers
  refseq.fasta.proteins <- refseq.id.length.df[,one.prot.colname]
  
  refseq.protein.ids.df[,one.prot.colname] <- refseq.fasta.proteins[pmatch(refseq.protein.ids.df[,"labkey_prot_name"],
                                                                           refseq.fasta.proteins)]
  
  refseq.id.length.index.df <- merge(refseq.id.length.df,
                                     refseq.protein.ids.df)
  
  prot.length.formula <- ""
  if(group.proteins.max.or.min == "max") {
    prot.length.formula <- as.formula(paste("~", "which.max(", prot.length.colname,")"))
  } else if(group.proteins.max.or.min == "min") {
    prot.length.formula <- as.formula(paste("~", "which.min(", prot.length.colname,")"))
  }
  
  refseq.id.protlength.df<- FilterColumnDplyrByGroup(input.df = refseq.id.length.index.df,
                                                     groupby.colnames = prot.colname,
                                                     input.formula = prot.length.formula)
  
  refseq.id.length.filt.df <- unique(refseq.id.protlength.df[,c(gene.colname,
                                                                prot.colname,
                                                                tx.colname,
                                                                prot.annot.colname,
                                                                one.prot.colname,
                                                                prot.length.colname)])
  return(refseq.id.length.filt.df)
}


# Normalizes for both protein length and total number of identified spectra from a sample
# NSAF = (SC/L)n/sum(SC/L)i
calcNSAFValues <- function(prot.df,
                           prot.exp.colname,
                           protseq.length.colname, 
                           nsaf.outfilename) {
  nsaf.calc.values <- prot.df[,prot.exp.colname]/prot.vars.results.sc.cbind[,protseq.length.colname]
  nsaf.values <- nsaf.calc.values/sum(nsaf.calc.values, na.rm = TRUE)
  
  return(nsaf.values)
}


processProteinFile <- function(){
  args <- commandArgs(TRUE)
  mzxml.file <- args[1]
  prot.file <- args[2]
  fasta.file <- args[3]
  output.file <- args[4]
  input.prot.colname <- args[5]
  sampleid.colname <- args[6]
  prot.exp.input.colname <- args[7]
  protseq.length.colname <- args[8]
  
  if(is.na(input.prot.colname)) {
    input.prot.colname <- "Protein"
  }
  
  if(is.na(sampleid.colname)) {
    sampleid.colname <- "sample_id"
  }

  if (is.na(prot.exp.input.colname)) {
    prot.exp.input.colname <- "Total.Filtered.Peptides"
  }
  
  if (is.na(protseq.length.colname)) {
    protseq.length.colname <- "prot_length"
  }
  
  output.prot.colname <- "proname"
  prot.exp.output.colname <- "SC"
  nsaf.colname <- "NSAF"
  
  
  # Brute force -> correct in future
  sample.id <- gsub("(.*)\\..*", "\\1", basename(mzxml.file))
  prot.df <- read.delim(prot.file, header = T, stringsAsFactors = FALSE)
  prot.df[,prot.exp.output.colname] <- prot.df[,prot.exp.input.colname]
  prot.df[,output.prot.colname] <- prot.df[,input.prot.colname]
  prot.df[,sampleid.colname] <- sample.id
  #Retain only required column names
  prot.filt.df <- unique(prot.df[,c(output.prot.colname,
                                    sampleid.colname,
                                    prot.exp.output.colname)])
  # Get protein lengths from Uniprot:
#   prot.filt.df <- getProteinLengthFromRefseqID(protein.fastafile = fasta.file,
#                                                refseq.protein.groups = prot.filt.df[,output.prot.colname],
#                                                prot.colname = output.prot.colname,
#                                                one.prot.colname = labkey.prot.colname,
#                                                prot.length.colname = protseq.length.colname,
#                                                group.proteins.max.or.min)
#   prot.filt.df[,nsaf.colname] <- calcNSAFValues(prot.filt.df,
#                                                 prot.exp.output.colname,
#                                                 protseq.length.colname, 
#                                                 nsaf.outfilename)
#   prot.nsaf.df <- unique(prot.df[,c(output.prot.colname, sampleid.colname,
#                                     output.prot.colname, nsaf.colname)])
#   print(prot.nsaf.df[1:5,])
  write.table(prot.filt.df, output.file, sep = "\t", row.names = F, quote = F)
  return(output.file)
}

processProteinFile()