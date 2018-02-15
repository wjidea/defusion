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

parser = argparse.ArgumentParser(description = 'process fasta files to drop fused gene entries')

# parser.add_argument('-i', '--gff_input', help='gff_file',required = True)
parser.add_argument('-i', '--fasta_in', help='fasta protein or transcript sequences',required = False)
parser.add_argument('-g', '--gff_in', help='gff file after defused genes',required = True)
parser.add_argument('-o', '--fasta_out', help='output fasta',required = True)

# parser.add_argument('-o', '--output', help = 'output AED scores', required = False)

parser.add_argument('-f', '--force', help = 'force overwrite', action = 'store_true')
parser.add_argument('-v', '--verbose', help = 'increase verbosity', action = 'store_true')


# extract mRNA name and aed scores from the gff file after defusion
# this gff file only contains the gene features survived after defusion step
def extract_geneid_AED(gff_file):
    """
    extract AED score from genome
    :param gff_file:
    :param goi_file:
    :return:
    """
    fi_gff = open(gff_file, 'r')
    aed_prog = re.compile("_AED=([^;]*)")
    gene_id_prog = re.compile("ID=(.*?)[;\n]")
    goi_list = []
    
    for line in fi_gff.readlines():
        line_L = line.rstrip().split('\t')
        
        if len(line_L) > 3 and line_L[2] == "mRNA":
            gene_id_search = gene_id_prog.search(line)
            aed_search = aed_prog.search(line)
            print gene_id_search.group(1), aed_search.group(1)
            goi_list.append(gene_id_search.group(1))
        else:
            continue
        
    fi_gff.close()
    
    return(goi_list)
    
    
def drop_fasta_entries(fasta_file, keep_list, fasta_out):
    seq_fn = SeqIO.parse(fasta_file, 'fasta')
    seq_fo = open(fasta_out, 'wb')
    for seqrec in seq_fn:
        if seqrec.id in keep_list:
            SeqIO.write(seqrec, seq_fo,'fasta')
    
    seq_fo.close()


def main():
    args = parser.parse_args()
    
    # check if the input fild existed
    input_gff = args.gff_in
    input_fasta = args.fasta_in
    out_file = args.fasta_out

    # goi_list = extract_AED_scores(gff_file=input_gff, goi_file=input_goi)
    goi_list = extract_geneid_AED(input_gff)
    drop_fasta_entries(input_fasta, goi_list, out_file)
    
if __name__ == "__main__":
    main()