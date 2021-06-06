#!/usr/bin/env python

# this script to sync gff features with transcript/protein sequences

import gffutils
import argparse
from Bio import SeqIO
from itertools import izip


def get_args():
    parser = argparse.ArgumentParser(description = 'Sync gff with fasta transcript or protein input')

    # parser.add_argument('-i', '--gff_input', help='gff_file',required = True)
    parser.add_argument('-t', '--trans_in', help='fasta transcript sequences',required = False)
    parser.add_argument('-p', '--prot_in', help='fasta protein sequences',required = False)
    parser.add_argument('-g', '--gff_in', help='gff.db sqlite file after defused genes',required = True)
    parser.add_argument('-o', '--gff_out', help='output fasta', required = True)

    parser.add_argument('-f', '--force', help = 'force overwrite', action = 'store_true')
    parser.add_argument('-v', '--verbose', help = 'increase verbosity', action = 'store_true')
    args = parser.parse_args()
    return args

def get_target_gene_names(fasta_file):
    '''
    get gene names list from input fasta file
    '''
    seq_fn = SeqIO.parse(fasta_file, 'fasta')
    gene_list = []
    for seqrec in seq_fn:
        gene_list.append(seqrec.name.rstrip('-mRNA-1'))
    print('Total genes in fasta: {}'.format(len(set(gene_list))))
    return gene_list

def rm_gff_features(gff_in, gff_out, gene_list):
    """
    rm gff features that are not in the fasta file
    Args:
        gff_in ([sqlite file]): [description]
        gff_out ([gff text]): [description]
        gene_list ([list]): [description]
    """
    db = gffutils.FeatureDB(gff_in)
    gff_fo = open(gff_out, 'w')
    gff_fo.write('## gff-version 3\n')
    for feature in db.all_features():
        if feature.featuretype == 'contig':
            gff_fo.write(str(feature) + '\n')
        # elif feature.source == 'maker':
        if feature.featuretype == 'gene' and feature.id in gene_list:
            gff_fo.write(str(feature) + '\n')
            for child_feature in db.children(feature.id, order_by='featuretype'):
                gff_fo.write(str(child_feature) + '\n')
        if feature.source == 'repeatmasker':
            gff_fo.write(str(feature) + '\n')
    gff_fo.close()

def validate_fasta_gff(gffdb, gene_list):
    """
    determine the gene list that intersect between gff and fasta files
    Args:
        gffdb ([type]): [description]
        gene_list ([type]): [description]

    Returns:
        [type]: [description]
    """
    gene_not_in_gff = []
    db = gffutils.FeatureDB(gffdb)
    for gene in gene_list:
        try:
            db[gene]
        except:
            gene_not_in_gff.append(gene)
    print('Number of genes that are not in gff file: {}'.format(len(gene_not_in_gff)))
    return gene_not_in_gff

def rm_fasta_records(rm_gene_list, prot_fa, trans_fa, suffix='-mRNA-1'):
    trans_fh = SeqIO.parse(trans_fa, 'fasta')
    prot_fh = SeqIO.parse(prot_fa, 'fasta')

    trans_out = []
    prot_out = []

    for trans, prot in izip(trans_fh, prot_fh):
        if trans.name.rstrip(suffix) not in rm_gene_list:
            trans_out.append(trans)
            prot_out.append(prot)

    SeqIO.write(trans_out, 'trans_synced.fasta', 'fasta')
    SeqIO.write(prot_out, 'protein_synced.fasta', 'fasta')

def main():
    args = get_args()
    gene_list = get_target_gene_names(args.trans_in)
    rm_gff_features(args.gff_in, args.gff_out, gene_list)
    rm_gene_list = validate_fasta_gff(args.gff_in, gene_list)
    rm_fasta_records(rm_gene_list, args.prot_in, args.trans_in)

if __name__ == '__main__':
    main()
