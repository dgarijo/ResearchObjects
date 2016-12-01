args <- commandArgs(TRUE)

library(customProDB)

print(args)
## Run the function to generate the protein sequence fasta file
easyRun(bamFile = args[1],
        vcfFile = args[2],
        outfile_path = dirname(args[3]),
        outfile_name = basename(args[3]),
        annotation_path = args[4],
        INDEL = args[5],        
        COSMIC = args[6],
        lablersid = args[7],
        bedFile = args[8])