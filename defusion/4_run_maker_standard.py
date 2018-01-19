#! /usr/bin/python

# run_maker_standard.py

# Jan. 18 2018
# Jie Wang
# Last update: Jan. 18, 2018

# description: run maker standard
# MAKER standard is to filter out the gene features that does not contain Pfam domain or a less than 1 AED score.

# input MAKER gff and fasta files
# output MAKER standard gff and fasta files

# this python script is a remake of Kevin Childs' maker standard pipeline in Perl

import re
import os, sys
import shlex
import argparse
import logging
import subprocess
from Bio import SeqIO
from utility_functions import input_validate, which, fix_path_slash


parser = argparse.ArgumentParser(description = 'run MAKER standard on gff and fasta files')

parser.add_argument('-t', '--trans_in', help='transcript sequence file [fasta]',required = True)
parser.add_argument('-p', '--prot_in', help='protein sequence file [fasta]',required = True)
parser.add_argument('-g', '--gff_in', help='gff file after MAKER annotation',required = True)
parser.add_argument('-o', '--output_dir', help='output directory',required = True)
parser.add_argument('-a', '--pfam', help='Pfam-A hmm file',required = True)
parser.add_argument('--pfam_cutoff', help='hmmer pfam search evalue cutoff',required = False, default="1e-10")

# parser.add_argument('-o', '--output', help = 'output AED scores', required = False)

parser.add_argument('-f', '--force', help = 'force overwrite', action = 'store_true')
parser.add_argument('-v', '--verbose', help = 'increase verbosity', action = 'store_true')

args = parser.parse_args()

tran_in = args.trans_in
prot_in = args.prot_in
gff_in = args.gff_in
pfam_in = args.pfam
out_dir = fix_path_slash(args.output_dir)
verbose = args.verbose
hmmer_cutoff = args.pfam_cutoff

if verbose:
    logging.basicConfig(level=logging.DEBUG, format='(%(processName)-10s) %(asctime)s %(message)s')
else:
    logging.basicConfig(level=logging.INFO, format='(%(processName)-10s) %(asctime)s %(message)s')

input_file_list = [tran_in, prot_in, gff_in]
input_validate(input_file_list)


# make tmp directory
try:
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    # if not os.path.exists(seqDir):
    #     os.mkdir(seqDir)
except OSError:
    print('Error in the temp directory')
    sys.exit()

# check HMMER executable in $PATH
if not which("hmmscan"):
    logging.error('hmmscan is not in path\n module load HMMER\n')
    sys.exit()
else:
    hmmer_scan = which("hmmscan")
    logging.info("HMMER is loaded")


def run_hmmer(prot_fasta, pfam):
    
    logging.info('Run HMMER on protein sequences')
    pfam_out = out_dir + "pfam_out.txt"
    
    if os.path.exists(pfam_out):
        return(pfam_out)
    else:
        cmd1 = '{0} --domE 1e-5 -E 1e-5 --cpu 4 --tblout {1} {2} {3}'.format(hmmer_scan, pfam_out, pfam, prot_fasta)
        cmd1L = shlex.split(cmd1)
        print cmd1
        p1 = subprocess.Popen(cmd1L, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p1.communicate()
        if stderr:
            logging.error(stderr)
        else:
            logging.info(stdout)
        return(pfam_out)
        
    
def generate_maker_standard_gene_list(hmmer_out, hmmer_cutoff, gff_in):
    
    pfam_std_set = set()
    maker_std_list = []

    fasta_prog = re.compile("#FASTA")
    aed_prog = re.compile("_AED=([^;]*)")
    gene_id_prog = re.compile("ID=(.*?)[;\n]")
    
    fh_hmmer = open(hmmer_out, 'rb')
    fh_gff = open(gff_in, 'rb')

    
    # parse hmmer output
    for line in fh_hmmer.readlines():
        
        if line.startswith("#"):
            continue
        
        # replace spaces with a tab
        line = re.sub("\s+", "\t", line)
        lineL = line.rstrip().split("\t")
        
        gene_id = lineL[2]
        evalue = float(lineL[4])
        
        if evalue <= float(hmmer_cutoff):
            pfam_std_set.add(gene_id)

    # parse gff to get aed score
    for line in fh_gff.readlines():
        if line.startswith("#"):
            continue
        if fasta_prog.search(line):
            print("match #fasta")
            return(maker_std_list)
        
        lineL = line.rstrip().split()
        
        if lineL[1] == 'maker' and lineL[2] == 'mRNA':
            feat_id_search = gene_id_prog.search(lineL[8])
            if feat_id_search:
                feat_id = feat_id_search.group(1)
                aed_search = aed_prog.search(lineL[8])
                if aed_search:
                    aed_score = float(aed_search.group(1))
                    # print(feat_id, aed_score)
                    if aed_score < 1.0 or feat_id in pfam_std_set:
                        maker_std_list.append(feat_id)
                else:
                    print("failed to parse aed score")
                    sys.exit()
            else:
                print("failed to parse gene id")
                sys.exit()

    fh_hmmer.close()
    fh_gff.close()

    logging.info("Maker standard gene: {}".format(len(maker_std_list)))
    return(maker_std_list)


def convert_gene_id(gene_list):
    # convert mRNA ID to gene ID in gff file
    gene_id_list = []
    for mRNA_name in gene_list:
        gene_id = re.sub("-mRNA-\d+$", "", mRNA_name)
        gene_id_list.append(gene_id)
    return(gene_id_list)

def create_maker_standard_gff(gene_list, gff):
    
    logging.info("Start gff standard process")
    
    mRNA_name_set = set()
    
    # compile regex
    out_gff_name = out_dir + re.sub("gff$", "std.gff", gff)
    
    fasta_prog = re.compile("#FASTA")
    gene_name_prog = re.compile("Name=(.*?)[;\n]")
    paranet_name_prog = re.compile("Parent=(.*?);Name=(.*?)[;\n]")
    parent_id_prog = re.compile("Parent=(.*?)[;\n]")
    
    fh_gff = open(gff, "rb")
    fo_gff = open(out_gff_name, "wb")
    
    for line in fh_gff:
        if line.startswith("#"):
            fo_gff.write(line)
            continue
            
        if fasta_prog.search(line):
            print("match #fasta")
            return(None)
        
        lineL = line.rstrip().split()
        
        # keep all non-maker features
        if lineL[1] != "maker":
            fo_gff.write(line)
            continue
        
        # keep gene feature
        if lineL[2] == "gene":
            gene_name_search = gene_name_prog.search(line)
            if gene_name_search:
                gene_name = gene_name_search.group(1)
                if gene_name in gene_list:
                    fo_gff.write(line)
                    continue
            
        # keep mRNA feature
        if lineL[2] == "mRNA":
            paranet_name_search = paranet_name_prog.search(line)
            if paranet_name_search:
                parent_id = paranet_name_search.group(1)
                mRNA_name = paranet_name_search.group(2)
                if parent_id in gene_list:
                    fo_gff.write(line)
                    mRNA_name_set.add(mRNA_name)
                    continue
        
        # keep exon CDS UTR features
        if lineL[2] == "exon" or lineL[2] == "CDS" or lineL[2] == "UTR":
            paranet_search = parent_id_prog.search(line)
            if paranet_search:
                parent_id = paranet_search.group(1)
                if parent_id in mRNA_name_set:
                    fo_gff.write(line)
                    continue

    fo_gff.close()
    fh_gff.close()

    logging.info("Finish gff standard process")
    
def create_maker_standard_fasta(gene_list, fasta):
    logging.info("Start fasta standard process")
    seq_fn = SeqIO.parse(fasta, 'fasta')
    prefix, ext = os.path.splitext(fasta)
    out_fasta_file = out_dir + prefix + ".std.fa"
    seq_fo = open(out_fasta_file, 'wb')
    
    for seqrec in seq_fn:
        if seqrec.id in gene_list:
            SeqIO.write(seqrec, seq_fo, 'fasta')
    seq_fo.close()
    logging.info("Finish fasta standard process")


def main():
    pfam_out = run_hmmer(prot_in, pfam_in)
    std_gene_list = generate_maker_standard_gene_list(pfam_out, hmmer_cutoff, gff_in)
    gene_id_list = convert_gene_id(std_gene_list)
    
    create_maker_standard_gff(gene_id_list, gff_in)
    create_maker_standard_fasta(std_gene_list, prot_in)
    create_maker_standard_fasta(std_gene_list, tran_in)
    
    
if __name__ == "__main__":
    main()




