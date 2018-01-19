#! /usr/bin/python

# process_fasta

# Jan 08, 2018
# Jie Wang
# Last update: Jan 08, 2018


# remove the split genes from the old protein fasta files using brk format


import re

from argparse import ArgumentParser
from Bio import SeqIO

####### Header file ########
description_arg = ''
usage_arg = 'python process_fasta.py'

parser = ArgumentParser(description=description_arg, usage=usage_arg)

parser.add_argument('-i', '--seq', help='sequence in fasta format', required=True)
parser.add_argument('-b', '--brk', help='brk file', required=True)
parser.add_argument('-o', '--out', help='output', required=True)
args = parser.parse_args()

in_seq = args.seq
out_seq = args.out
in_brk = args.brk

seq_fh = SeqIO.index(in_seq, 'fasta')

def parse_brk(brk):
    geneList = []
    with open(brk, 'rb') as fh:
        for line in fh.readlines():
            lineL = line.rstrip().split()
            gene_id = lineL[1]
            geneList.append(gene_id)
    return(geneList)

gene_list = parse_brk(in_brk)

with open(out_seq, 'wb') as fo:
    for seq_id in seq_fh:
        if seq_id not in gene_list:
            # print(seq_id, str(seq_fh[seq_id].seq))
            fo.write(">{}\n{}\n".format(seq_fh[seq_id].id, str(seq_fh[seq_id].seq)))
            


        