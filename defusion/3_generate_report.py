#! /usr/bin/python

# generate_report.py

# Aug 30, 2017
# Jie Wang
# Last update: Jan. 10, 2018


# description: extract the AED score from the original and improved transcript sequence header
# plot them in matplotlib and seaborn

import re
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description = 'parse BLAST output and write to tab-delimited file')

parser.add_argument('-i', '--gff_input', help='gff_file',required = True)
parser.add_argument('-s', '--fasta_input', help='fasta protein or transcript sequences',required = False)
parser.add_argument('-g', '--goi_input', help='breaking coordinate file',required = True)
parser.add_argument('-o', '--output', help = 'output AED scores', required = False)

parser.add_argument('-f', '--force', help = 'force overwrite', action = 'store_true')
parser.add_argument('-v', '--verbose', help = 'increase verbosity', action = 'store_true')


def extract_fasta_seqids(fasta_file):
    ini = 0


def extract_AED_scores(gff_file, goi_list):
    """
    extract AED score from genome
    :param gff_file:
    :param goi_file:
    :return:
    """
    fi_gff = open(gff_file, 'r')
    aed_prog = re.compile("_AED=([^;]*)")
    gene_id = re.compile("ID=([^;]*);")
    
    for line in fi_gff.readlines():
        line_L = line.rstrip().split('\t')
        if len(line_L) > 3 and line_L[2] == "mRNA":
        
            attr = line_L[8]
            gene_id = gene_id.search()
        
            if gene_id in goi_list:
                print line,
        else:
            continue
        
    fi_gff.close()
    
    return(goi_list)
    
    
def drop_fasta_entries(fasta_file, drop_list):
    seq_fn = SeqIO.parse(fasta_file, 'fasta')
    seq_fo = open("seq_out.fasta", 'wb')
    for seqrec in seq_fn:
        if seqrec.id not in drop_list:
            SeqIO.write(seqrec, seq_fo,'fasta')
    
    seq_fo.close()


def main():
    args = parser.parse_args()
    
    # check if the input fild existed
    input_gff = args.gff_input
    input_goi = args.goi_input
    input_fasta = args.fasta_input
    out_file = args.output

    goi_list = extract_AED_scores(gff_file=input_gff, goi_file=input_goi)
    drop_fasta_entries(input_fasta, goi_list)
    
if __name__ == "__main__":
    main()
    