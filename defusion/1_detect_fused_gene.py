#! /usr/bin/python

# defuse_gene.py

# Oct. 14 2016
# Last update: Feb. 15, 2018
# Jie Wang

# description: Identify the mis-annotated genes in the MAKER annotation pipeline,
# particularly for the tandem duplicated genes annotated as fused gene in the
# annotation project.
#
# usage:
#     input1: -i input transcript nucleotide sequence fasta file
#     input2:  ?? configuration file ??
#     output1: list of fused gene identifier
#     outpur2: estimate the breaking point for the duplicated gene
#
# implementation steps:
# 1, read in the transcript DNA sequence and seprate them into individual file
# 2, run self BLAST in parallel
# 3, write report in blastn out format 6
#
# This output will provide the out for potential fused gene candidates


import os
import sys
import glob
import logging
import multiprocessing
import gffutils

from utility_functions import which
from argparse import ArgumentParser
from shutil import copyfile
from Bio import SeqIO
from functools import partial # used in multipleprocess manager lock, pass lock to function
from Bio.Blast.Applications import NcbiblastnCommandline


####### Header file ########

description_arg = 'First step to detect the tandem duplicated fused gene annoatation'
usage_arg = 'python 1_detect_fused_gene -h to show help information'

parser = ArgumentParser(description=description_arg, usage=usage_arg)
parser.add_argument('-i', '--input', help='transcripts in [fasta] format', required=True)
parser.add_argument('-g', '--gffin', help='input gff from MAKER2 output', required=True)
parser.add_argument('-d', '--gffdb', help='input gffdb sqlite from previous run', required=False, default=None)
parser.add_argument('-n', '--numOfProcess', help='number of processes to be included', default=4, required=False)
parser.add_argument('-b', '--blastn', help='path to executable blastn', default='', required=False)
parser.add_argument('-p', '--prefix', help='prefix dirctory', default='tmp/', required=False)
parser.add_argument('-f', '--force', help='force overwrite existing output', action='store_true', required=False)
parser.add_argument('-v', '--verbose', help='increase verbosity', action='store_true')
args = parser.parse_args()

# check if the input fild existed
inFile = args.input
# outPath = args.outpath
# outBlast = args.blastOut
numOfProcess = int(args.numOfProcess)
blastn = args.blastn
gffIn = args.gffin
gffdb = args.gffdb
prefixDir = os.path.join(args.prefix,'')
seqDir = os.path.join(prefixDir+'seqs', '')
inFilesL = [inFile]
# outFilesL = []

# TODO add a function to validate input file, such as fasta and gff

if args.verbose:
    logging.basicConfig(level=logging.DEBUG, format='(%(processName)-10s) %(asctime)s %(message)s')
else:
    logging.basicConfig(level=logging.INFO, format='(%(processName)-10s) %(asctime)s %(message)s')

if not os.path.isfile(gffIn):
    logging.error('gff file is not in the specified path')
    sys.exit()

# check executables in $PATH

if not blastn: # test makerbin is provided
    if not which("blastn"):
        logging.error('blastn executable is not in path\n module load BLAST+\n')
        sys.exit()
    else:
        blastn = which("blastn")
        logging.info("blastn is loaded in $PATH")


# check input file list
for inputFile in inFilesL:
    if not os.path.isfile(inputFile):
        logging.error('Input file {0} does not exist!'.format(inputFile))
        sys.exit()

# # check if output file is present
if args.force:
    pass
else:
    if os.path.isdir(prefixDir):
        print('output file {0} already exists'.format(prefixDir))
        sys.exit()

# make tmp directory
try:
    if not os.path.exists(prefixDir):
        os.mkdir(prefixDir)
    if not os.path.exists(seqDir):
        os.mkdir(seqDir)
except OSError:
    logging.error('Error in the temp directory')
    sys.exit()

######## scripts start from here #########

def seprate_fasta(inFile):
    # outPath = outPath.rstrip('/')
    outPath = prefixDir + 'seqs/'
    allSeq = SeqIO.parse(inFile, 'fasta')
    for seqRec in allSeq:
        SeqIO.write(seqRec, '{0}{1}.fasta'.format(outPath, seqRec.id), 'fasta')


def self_BLAST(lock, seq):
    """
     Run self BLASTn and report seqID if fits fused gene annotation patterns
    :param lock: process lock for running BLAST in parallel
    :param seq: sequence file handle
    :return: sequence ID that meet the requirement
    """
    blastfmt = "'6 qseqid qlen length sstart send qstart qend evalue'"
    blastnCmd = NcbiblastnCommandline(cmd=blastn, query=seq, subject=seq,
                                      task='blastn', evalue='1e-5',
                                      outfmt=blastfmt, max_hsps='10')
    blastnResults = blastnCmd()
    blastHits = blastnResults[0].rstrip().split('\n')
    blastHitsSplit = [hit.split('\t') for hit in blastHits]
    
    # more than 2 hits and aln length large than a 1/4 of seq length
    logging.debug('start running selfblast, {}'.format(seq))
    
    if len(blastHits) > 2:
        if float(blastHitsSplit[1][2]) > (float(blastHitsSplit[0][1]) / 4.0):
            with lock:
                with open(prefixDir + 'selfBlastOut.txt', 'ab') as handle:
                    
                    if len(blastHits) < 4:
                        handle.write(blastHitsSplit[0][0] + '\t2\n')
                    elif len(blastHits) > 3:
                        handle.write(blastHitsSplit[0][0] + '\t3+\n')

                return(blastHitsSplit[0][0])
    else:
        return(None)
    
def run_blast_parallel(seqFileL):
    # instantiate jobs using mp.Pool
    pool = multiprocessing.Pool(numOfProcess)
    manager = multiprocessing.Manager()
    lock = manager.Lock()
    
    geneListFile=  prefixDir+'selfBlastOut.txt'
    fusedGeneL = []
    geneOutL = []
    
    print seqFileL
    
    if os.path.exists(geneListFile):
        # if selfBLAST step is done, read the results then
        with open(geneListFile, 'rb') as fh:
            for line in fh.readlines():
                lineL = line.rstrip().split()
                geneOutL.append(lineL[0])
    else:
        try:
            fusedGeneL = pool.map(partial(self_BLAST, lock), seqFileL)
            # fusedGeneL = fusedGeneL + [result for result in blastResults]
            pool.close()
            pool.join()
    
        except KeyboardInterrupt:
            # Allow ^C to interrupt from any thread.
            sys.stdout.write('\033[0m')
            sys.stdout.write('User Interupt\n')
        
        geneOutL = [gene_id for gene_id in fusedGeneL if gene_id is not None]
        
        logging.info('{}\n'.format(geneOutL))
    return(geneOutL)


def rm_gff_features(gff_in, prefixDir):
    '''
    rm features from exonerate or gene prediction tools
    Only keep maker
    :param gffin: input gff file
    :return: modified gff file without evidence
    '''
    logging.info("remove ununsed features from gff file to increase load speed")
    
    fh_gff = open(gff_in, 'rb')
    fo_gff = open(prefixDir + "modified_annotation.gff", 'wb')
    
    for line in fh_gff.readlines():
        
        if line.startswith("##FASTA"):
            break
        
        if line.startswith("##gff-version 3"):
            fo_gff.write(line)
            continue
        elif line.startswith("##"):
            continue
        
        lineL = line.rstrip().split()
        
        entry_source = lineL[1]
        
        if entry_source == "repeatmasker" or entry_source == "maker" or entry_source == r".":
            fo_gff.write(line)
    
    fo_gff.close()
    fh_gff.close()
    
    return (prefixDir + "modified_annotation.gff")
    

def filter_gff(geneL, gff):
    """
    :param geneL: list of candidate genes
    :param gff: maker gff out file
    :return: list of geneID that can need to be filtered
    """
    geneDict = {}
    dbName = prefixDir + 'gff.db'
    
    if gffdb and (not os.path.isfile(dbName)):
        copyfile(gffdb, dbName)
        
    try: # read gff to sqlite db
        if os.path.isfile(dbName):
            db = gffutils.FeatureDB(dbName)
            logging.debug('Found gff.db, now loading')
        else:
            logging.info('Loading gff to sqlite database, please wait.')
            db = gffutils.create_db(gff, dbfn=dbName, force=False, keep_order=False, merge_strategy='create_unique',
                                    id_spec=['ID'], sort_attribute_values=False)
            # create index for speed improvement, the original index within gffutils seems failed somehow
            db.execute('CREATE INDEX coordinate_index ON features(seqid, start, end)')
            logging.debug('complete loading gff.db')
    except ValueError:
        logging.error('import gff/gtf is not present in given path')
    
    for gene in geneL: # loop through the list and collect contig info
        gene_fused = db[gene]
        g_start = gene_fused.start
        g_end = gene_fused.end
        g_seqid = gene_fused.seqid
        
        # count how many exon in this gene
        g_exon_counter = 0
        logging.debug('count how many exons')
        for i in db.children(gene_fused, featuretype='exon'):
            g_exon_counter += 1
        
        if g_exon_counter <= 1:
            continue
        
        geneDict[gene] = [g_seqid, g_start, g_end, g_exon_counter] # seqid, start, end
    
    #fixme 1 if the gene lenth (including introns) is way longer than the blastn query hit? should I delete?
    #fixme 2 example: A07:23319901..23327780 (7.88 Kb); match length and gene length

    logging.debug('start writing gene coordinate file')
    with open(prefixDir+'geneCoordinate.txt', 'wb') as f:
        for key in geneDict:
            f.write('{0}\t{1}:{2}..{3}\t{4}\n'.format(key, *geneDict[key]))

    return(geneDict)

# def breakPont(geneDict, seqIn, blastOut):
def break_point(geneDict):
    """
    Parse blast fmt 6 output and extract parts of the sequence for re-annotation
    
    ---|------region of interest 1------|------region of interest 2------|---
    break point 1                    break point 2                  break point 3
    
    break point 1 = left border - 100
    break point 3 = right border + 100
    break point 2 = mid point between two matched features or mid point of gene
    
    Another big problem is to have > 2 tandem duplicates, which can complicate the problem
    
    Possible solution is to check the output transcripts to run the detect fused gene code recursively, until no fused
    annotation found in the transcripts file.
    
    For the un-determined the one, raise flag for manual inspections.
    
    :param geneDict: geneID:[seqid, start, end]
    :param seqIn: sequence file in fasta format
    :param blastOut: txt file generate from runBlastParallel
    :return: chunks of individual fasta sequence files
    """
    # determine the breaking point TODO improve breaking point selection composite score
    # with open(blastOut)
    
    # implement the simple version of this function divide them in the mid point
    dbName = prefixDir + 'gff.db'
    
    try:
        db = gffutils.FeatureDB(dbName)
        logging.debug('Found gff.db, now loading')
        # dbe = sqlite3.connect(dbName)
    except ValueError:
        print('Check gff.db present in working directory')
        sys.exit()
        
    blast_dict = {}
    with open(prefixDir + 'selfBlastOut.txt', 'rb') as handle:
        for line in handle.readlines():
            lineL = line.rstrip().split()
            blast_dict[lineL[0]] = lineL[1]
    
    with open(prefixDir+'break_coordinates.brk', 'wb') as fh:
        for gene in geneDict:
            geneFeatL = geneDict[gene]
            # ROI = db.region(seqid=geneFeatL[0], start=geneFeatL[1], end=geneFeatL[2], completely_within=False)
            seqID = geneFeatL[0]
            seqStart = geneFeatL[1] - 50 # increase the boundaries to expand the two borders
            seqEnd = geneFeatL[2] + 50
            # this is simple version just to demonstrate the format settings
            point2Break = int((int(seqEnd) - int(seqStart))/2) + int(seqStart)
            lineStr = '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(seqID, gene, seqStart, seqEnd, point2Break, blast_dict[gene])
            fh.write(lineStr)
        # todo more accurate version for determining mid point need to learn more the BLAST results to detemine how
        # many repeats in the regions
        
def chimera_fuse():
    """
    chimera fused gene annotation in MAKER, especially with SringTie Associated genes
    http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-279
    
    KOG database or just to use orthoDB database, another option
    :return:
    """
    init = 0


def main():
    seprate_fasta(inFile)
    seqFileL = glob.glob(prefixDir + 'seqs/*.fasta')

    geneL = run_blast_parallel(seqFileL)
    
    modGFF = rm_gff_features(gffIn, prefixDir)
    geneDict = filter_gff(geneL, modGFF)
    break_point(geneDict)

if __name__ == '__main__':
    main()