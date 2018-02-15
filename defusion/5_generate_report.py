#! /usr/bin/python

# generate aed improvement score

# Jan 08, 2018
# Jie Wang
# Last update: Feb 14, 2018


# extract AED scores from old and defused gffs files


import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from argparse import ArgumentParser

####### Header file ########
description_arg = ''
usage_arg = 'python 5_generate_report -i original_MAKER.gff -g merged_defused.all.mod.gff ' \
            '-b curated_breaking_point.brk -o output_prefix'

parser = ArgumentParser(description=description_arg, usage=usage_arg)

parser.add_argument('-i', '--gold', help='Original MAKER gff file', required=True)
parser.add_argument('-g', '--gdef', help='defused GFF file', required=True)
parser.add_argument('-b', '--brk', help='brk file', required=True)
parser.add_argument('-o', '--out', help='output prefix for output files', required=True)
args = parser.parse_args()


# extract mRNA name and aed scores from the gff file after defusion
# this gff file only contains the gene features survived after defusion step
def extract_geneid_AED(gff_file, sublist = None):
    """
    extract AED score from genome
    :param gff_file:
    :param goi_file:
    :return:
    """
    fi_gff = open(gff_file, 'r')
    aed_prog = re.compile("_AED=([^;]*)")
    gene_id_prog = re.compile("ID=(.*?)[;\n]")
    aed_list = []
    
    for line in fi_gff.readlines():
        line_L = line.rstrip().split('\t')
        
        if len(line_L) > 3 and line_L[2] == "mRNA":
            gene_id_search = gene_id_prog.search(line)
            aed_search = aed_prog.search(line)
            # print gene_id_search.group(1), aed_search.group(1)
            gene_id = gene_id_search.group(1)
            aed_score = aed_search.group(1)
            # goi_list.append(gene_id_search.group(1))
            
            if sublist:
                if gene_id in sublist:
                    aed_list.append(aed_score)
                else:
                    continue
            else:
                aed_list.append(aed_score)
    fi_gff.close()
    aed_list.sort()
    
    return(aed_list)


def parse_brk(brk_file):
    fh = open(brk_file, 'rb')
    
    gene_id_list = []
    for line in fh.readlines():
        lineL = line.rstrip().split()
        gene_id = lineL[1]
        gene_id_list.append(gene_id)
    fh.close()
    
    return(gene_id_list)


def bin_aed_score(sorted_aed_list, outfile, label, num_bin = 40):
    
    sorted_aed_list = map(float, sorted_aed_list)
    hist, edges = np.histogram(sorted_aed_list, bins=num_bin, range=[0,1])
    edges = np.delete(edges, -2)
    
    col1 = edges
    col2 = np.cumsum(hist)
    col3 = np.cumsum(hist) / float(sum(hist))
    
    with open(outfile, 'wb') as fh:
        for i in range(num_bin):
            print(col1[i], col2[i], col3[i])
            out_str = "{:}\t{}\t{}\t{}\n".format(col1[i], col2[i], col3[i], label)
            fh.write(out_str)


def draw_aed_plot(fused_file, defused_file, out_prefix):
    
    fused_aed_data = pd.read_table(fused_file, header=None, names=['aed', 'gene_counts', 'norm_gene', 'category'])
    defused_aed_data = pd.read_table(defused_file, header=None, names=['aed', 'gene_counts', 'norm_gene', 'category'])

    fig, ax = plt.subplots()

    color_fuse, color_defuse_std = 'blue', 'red'

    line1_fused, = ax.plot(fused_aed_data.aed, fused_aed_data.norm_gene, color=color_fuse)
    line2_defuse_std, = ax.plot(defused_aed_data.aed, defused_aed_data.norm_gene, color=color_defuse_std)

    plt.legend([line1_fused, line2_defuse_std], ['fused_genes', 'defused_genes'], loc=4, fontsize=14)

    ax.set(xlim=[0, 1.1])
    ax.set_xlabel('Annotation edit distance (AED)')
    ax.set_ylabel('Percentage of genes [0-1]')
    ax.set_title('AED scores improvement after defusion - {}'.format(out_prefix), fontweight="bold",
                 fontsize=14, y=1.0)

    # plt.show()
    plt.savefig('{}_AED_score.png'.format(out_prefix), dpi=600)
    
    
def main():
    old_gff = args.gold
    defuse_gff = args.gdef
    in_brk = args.brk
    prefix = args.out
    
    fused_gene_id = parse_brk(in_brk)
    fused_aed = extract_geneid_AED(old_gff, fused_gene_id)
    defused_aed = extract_geneid_AED(defuse_gff)

    fused_aed_file = "{}_fused_aed.txt".format(prefix)
    defused_aed_file = "{}_defused_aed.txt".format(prefix)
    
    bin_aed_score(defused_aed, outfile=defused_aed_file, label="defused")
    bin_aed_score(fused_aed, outfile=fused_aed_file, label="fused")

    draw_aed_plot(fused_aed_file, defused_aed_file, prefix)
    
if __name__ == "__main__":
    main()