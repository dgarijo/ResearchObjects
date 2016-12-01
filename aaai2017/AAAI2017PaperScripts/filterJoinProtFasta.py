### Compares two fasta files and generates a new joint fasta file with all unique proteins
### In case of overlaps, protein sequences from file1 are retained
import argparse
import re
from Bio import SeqIO

def generateFiltFastaFile(fasta1, fasta2, output_fasta):
    fasta1_records = SeqIO.parse(fasta1, "fasta")
    # Grab only the refseq ID from the fasta1-header in case we have additional id info
    fasta1_records_list = [one_fasta1_record for one_fasta1_record in fasta1_records]
    regex_refseqid = re.compile('(NP_\d+)')
    fasta1_protids_list = [prot_id.group(1) for header_id \
                            in fasta1_records_list \
                            for prot_id in [regex_refseqid.search(header_id.id)] \
                            if prot_id]

    fasta2_records = SeqIO.parse(fasta2, "fasta")
    # Remove any records from file2 (mostly reference proteome) that overlap with file1
    fasta2filt_records_list = [one_fasta2_record for one_fasta2_record \
                                in fasta2_records \
                                if one_fasta2_record.id not in fasta1_protids_list]
    fasta1_fasta2filt_records_list = fasta1_records_list + fasta2filt_records_list

    # Write the resulting fasta file
    output_fasta_handle = open(output_fasta, "w")
    SeqIO.write(fasta1_fasta2filt_records_list, output_fasta_handle, "fasta")
    output_fasta_handle.close()
    return output_fasta

def main():
    fasta1, fasta2, output_fasta = readAndParseCommandlineArgs()
    output_fasta = generateFiltFastaFile(fasta1, fasta2, output_fasta)
    print "Output file: ", output_fasta

#
# Read command line options
#
def readAndParseCommandlineArgs():
    parser = argparse.ArgumentParser(add_help=True, usage='%(prog)s [options]')
    parser.add_argument("-1", "--fasta1", \
                        dest="fasta1", action="store", default=False)
    parser.add_argument("-2", "--fasta2", \
                        dest="fasta2", action="store", default=False)
    parser.add_argument("-o", "--output", dest='out_fasta',\
                        action="store",\
                        default="filt_output.fasta")

    try:
        args = parser.parse_args()
        #print args
        return args.fasta1, args.fasta2, args.out_fasta
    except IOError, msg:
        parser.error(str(msg))

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()
