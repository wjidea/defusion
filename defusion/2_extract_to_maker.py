#! /usr/bin/python

# extract_to_maker.py

# Nov. 28 2016
# Jie Wang
# Last update: Jan. 18, 2018

# description: extract the sequence from the reference genome fasta file based on
# the gene coordinates generated from 01_detect fused gene step, and run MAKER
# final on each individual pieces of sequence. Use the Maker Final gff file to
# update the existing maker genome annotation files

# implementation steps:
# 1, extract the fasta based on the fused gene feature coordinates
# 2, write a wrap for running MAKER with all the provided pathes to the

import logging
import sys
import os
import subprocess
import shlex
import fileinput
import glob
import gffutils
import multiprocessing


from shutil import copy2
from argparse import ArgumentParser
from functools import partial
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utility_functions import which, fix_path_slash, log_err, parse_SeqID, \
    bioperl_loaded, check_coord, coord_error

####### Header file ########
description_arg = 'This is the second step of the defusion pipeline: 1) remove the original fused annotation' \
                  '2) run new MAKER session reannotate the fused annotation region'
usage_arg = 'python 2_extract_to_maker.py -h'

parser = ArgumentParser(description=description_arg, usage=usage_arg)

parser.add_argument('-i', '--genome', help='genome in fasta format', required=True)
parser.add_argument('-d', '--gffdb', help='input gff from MAKER2 output ', required=True)
parser.add_argument('-c', '--coordin', help='gene coodinate at break point', required=True)
parser.add_argument('-t', '--control', help='path trol file', required=True)
parser.add_argument('-m', '--makerbin', help='path to maker executables directory', required=False)
parser.add_argument('-n', '--numprocess', help='number of processes to be included', default=2, required=False)
parser.add_argument('-p', '--prefix', help='prefix dirctory', default='tmp/', required=False)
parser.add_argument('-s', '--datastore', help='with this flag to keep the datastore folder but need more disk space',
                    action='store_true', required=False, default=False)
parser.add_argument('-v', '--verbose', help='increase verbosity', action='store_true', default=False)
args = parser.parse_args()

seqIn = os.path.abspath(args.genome)
gffdb = os.path.abspath(args.gffdb)
ctlPath = fix_path_slash(args.control)
makerBinPath = args.makerbin
coordIn = args.coordin # breakfile
numProcess = int(args.numprocess)
prefixDir = fix_path_slash(args.prefix)
inFilesL = [seqIn, gffdb]
keep_datastore = args.datastore
verbose = args.verbose


if verbose:
    logging.basicConfig(level=logging.DEBUG, format='(%(processName)-10s) %(asctime)s %(message)s')
else:
    logging.basicConfig(level=logging.INFO, format='(%(processName)-10s) %(asctime)s %(message)s')

# make tmp directory
try:
    if not os.path.exists(prefixDir):
        os.mkdir(prefixDir)
    # if not os.path.exists(seqDir):
    #     os.mkdir(seqDir)
except OSError:
    print('Error in the temp directory')
    sys.exit()


# check input files
for inputFile in inFilesL:
    if not os.path.isfile(inputFile):
        print('Input file {0} does not exist!'.format(inputFile))
        sys.exit()

# check executables in $PATH

if not makerBinPath: # test makerbin is provided
    if not which("maker"):
        logging.error('MAKER executable is not in path\n module load MAKER\n')
        sys.exit()
    else:
        maker_exe_path = which("maker")
        makerBinPath, fname = os.path.split(maker_exe_path)
        makerBinPath = fix_path_slash(makerBinPath)
        logging.info("MAKER is loaded in $PATH {}".format(makerBinPath))

if bioperl_loaded:
    logging.info("BioPerl is loaded")
else:
    logging.error("BioPerl is not loaded")
    sys.exit()

######## scripts start from here #########
def flat_multiple_breaks(breakfile, outfile):
    fo = open(prefixDir + outfile, 'wb')
    coord_error_bool = False
    with open(breakfile, 'r') as fn:
        for line in fn.readlines():
            lineL = line.rstrip().split('\t')
            
            # set exception in the brk file
            lineL = [i for i in lineL if i] # rm empty item
            if lineL[2] == "collapsed" or len(lineL[2:]) < 3:
                continue # skip to next index
                
            coords_list = lineL[2:]
            coords_list = map(int, coords_list)
            coords_list.sort()
            
            for idx in range((len(coords_list) - 1)):
                coord_start = coords_list[idx]
                coord_end = coords_list[idx + 1]
                out_list = lineL[0:2] + [coord_start, coord_end]
                out_list = map(str, out_list)
                out_line = '\t'.join(out_list) + '\n'
                fo.write(out_line)
                
                # check if gene range large than 100k
                error_in_coord = check_coord(prefixDir, lineL[0], coord_start, coord_end)
                if not error_in_coord:
                    coord_error_bool = True
    
    # validate error presence and decide ignore or not
    coord_error(coord_error_bool)
    
    fo.close()


def extract_seq(line):
    """
    read input breakpoint file and extract sequence from the genome file
    :param line: line read from breakpoint file
    :return: a list to two folders
    """
    genomeSeq = SeqIO.index(seqIn, 'fasta')
    
    lineL = line.rstrip().split()
    seqID, geneID = lineL[0:2]
    start, end = lineL[2:4]
    folderL = [seqID, start, end]
    folder = '_'.join(folderL)
    
    start, end = map(int, lineL[2:4])
    seqSelect = genomeSeq[seqID][start:end].seq
    folderPath = os.path.join(prefixDir + folder, '')
    seqSelectFile = folderPath + folder + '.fa'
    
    try:
        logging.info('Create {}'.format(folderPath))
        os.mkdir(folderPath)
    except OSError:
        print('maker directory existed, please check')
        raise
    
    seq = open(seqSelectFile, 'wb')
    
    logging.info('Write sequence files {}'.format(seqSelectFile))
    seqrec = SeqRecord(seqSelect, folder, '', '')
    SeqIO.write(seqrec, seq, 'fasta')
    
    seq.close()
    
    return(folder)


def run_extract_seq_parallel(fileIn):
    
    with open(fileIn, 'rb') as fh:
        lineL = fh.readlines()
    
    folderL = []
    
    pool = multiprocessing.Pool(numProcess)
    
    try:
        folder = pool.map(extract_seq, lineL)
        folderL.extend(folder)
        pool.close()
        pool.join()
                
    except KeyboardInterrupt:
        # Allow ^C to interrupt from any thread.
        sys.stdout.write('\033[0m')
        sys.stdout.write('User Interupt\n')
    
    return(folderL)
    
    
def prep_maker_gff(gffdb, seqID, prefix):
    
    # os.chdir(prefix + seqID)
    # example seqID name Chr1_2741699_2750000
    logging.info('Start working on {} prep maker_gff'.format(seqID))
    
    # prepare the coordinates
    seqID_list = seqID.rstrip().split('_')
    
    end_ori = seqID_list.pop()
    start_ori = seqID_list.pop()
    chr = '_'.join(seqID_list)
    # chr, start_ori, end_ori = seqID.rstrip().split('_')
    start_ori = int(start_ori) - 1  # zero based
    end_ori = int(end_ori)
    end_new = end_ori - start_ori
    
    # load gffdb sqlite3
    try:
        db = gffutils.FeatureDB(gffdb)
    except:
        print('database problem')
        raise
    
    ROI = db.region(seqid=chr, start=start_ori, end=end_ori, completely_within=False)
    
    def adj_coord(old_coord, start_ori):
        old_coord = map(int, old_coord)
        new_coord = [coord - start_ori for coord in old_coord]
        if new_coord[0] < 1:
            new_coord[0] = 1
        if new_coord[1] > end_new:
            new_coord[0] = end_new
        return (new_coord)
    
    with open(prefix + '{0}/{0}.gff'.format(seqID), 'wb') as fh:
        fh.write('##gff-version 3\n')
        for feat in ROI:
            if feat.source in ['repeatmasker']:
                coord = [feat.start, feat.end]
                newCoord = adj_coord(coord, start_ori)
                feat.start, feat.end = newCoord[:]
                fh.write('{}\n'.format(str(feat)))
    
    # delete anything present within the range.
    # ROI_within = db.region(seqid=chr, start=start_ori, end=end_ori, completely_within=True)
    #
    # if os.path.isfile(gffdb + '.bak'):
    #     db.delete(ROI_within, make_backup=False)
    # else:
    #     db.delete(ROI_within, make_backup=True)
    logging.info('Finish {} prep maker_gff'.format(seqID))
    return({seqID:start_ori})

# extract sequence coordinate from breakpoint to delete all the associated gff features


def del_gene_records(gffdb, fileIn):
    '''
    select region of deletion and delete the entries
    :param gffdb: gff sqlite database
    :param fileIn: breakfile from step 1
    :return: no return
    '''
    logging.info('Start working on remove gene features from the analysis')
    
    # load gffdb sqlite3
    try:
        db = gffutils.FeatureDB(gffdb)
    except:
        logging.warning('error with sqlite database connection')
        raise
    
    # parse breakpoint file by line
    with open(fileIn, 'rb') as fh:
        for line in fh.readlines():
            lineL = line.split()
            seqID, gene_id = lineL[0:2]
            coordL = lineL[2:]
            coordL.sort()
            start = int(coordL[0]) - 1
            end = int(coordL[-1])
            
            print gene_id
            
            gene_id_feature = db[gene_id]
            
            print(gene_id_feature)
            
            strand_gene_id = gene_id_feature.strand
            # query1 = db.execute('''
            #           SELECT * FROM features WHERE start >= {0} AND end =< {1}
            #           '''.format(start, end))
            # print query1.fetchall()
            
            # extract region to be collected
            ROI_within = db.region(seqid=seqID, start=start, end=end, completely_within=True)
            # convert a generator to a list
            ROI_list = list(ROI_within)
            # get gene/mRNA ID from the mRNA entry
            geneIDs = [x.attributes['ID'][0] for x in ROI_list if x.featuretype == "mRNA"]
            
            keep_idx = []  # create a list container to save the index to keep in the original gff file or gffdb
            
            for num, i in enumerate(ROI_list):
                if i.featuretype == "exon" or i.featuretype == "CDS":
                    if i.attributes["Parent"][0] not in geneIDs:
                        keep_idx.append(num)
                elif i.strand != strand_gene_id:
                    keep_idx.append(num)
            
            # make the iterator to remove from gffdb
            ROI_list_rm = [element for i, element in enumerate(ROI_list) if i not in keep_idx]
            
            # remove associated features from the range
            if os.path.isfile(gffdb + '.bak'):
                db.delete(ROI_list_rm, make_backup=False)
            else:
                db.delete(ROI_list_rm, make_backup=True)
    
    logging.info('Complete remove fused gene annotation')


def prep_opt_ctl(seqID, prefix, ctlPath):
    # TODO check maker_opts.ctl if the files in the control files exists
    logging.info('[START] Run prep maker opt control file')
    seqPath = os.path.join(prefix + seqID, '')

    ctlFilesL = glob.glob(ctlPath + '*.ctl')

    for file in ctlFilesL:
        fileFullPath = os.path.abspath(file)
        copy2(fileFullPath, seqPath)  # copy file in shell

    optCtlFile = seqPath + 'maker_opts.ctl'
    seqFile = seqPath + seqID + '.fa'
    gffFile = seqPath + seqID + '.gff'
    fileInput = fileinput.input(optCtlFile, inplace=True, backup='.bak')
    try:
        if os.path.isfile(optCtlFile + '.bak'):
            print('bak file exists exit(0)')
            sys.exit()

        for line in fileInput:
            if line.startswith('genome='):
                print('genome={0}'.format(seqFile)),
            elif line.startswith('maker_gff='):
                print('maker_gff={0}'.format(gffFile))
            else:
                print(line),
    except IOError:
        print('There was a IO error')
        raise
    finally:
        fileInput.close()
        logging.info('[FINISH] prep maker opt control file')
     
def run_maker(seqID, prefix, makerbin, keepStore):
    os.chdir(prefix + seqID)
    
    # Run maker
    logging.info('start MAKER {}'.format(seqID))
    makerExe = makerbin + 'maker'
    
    cmd1 = '{maker} -q -fix_nucleotide'.format(maker=makerExe)
    cmd1L = shlex.split(cmd1)
    p1 = subprocess.Popen(cmd1L, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p1.communicate()
    log_err('MAKER', seqID, stdout, stderr)
    
    # Generate GFF file
    logging.info('start gff3_merge {}'.format(seqID))
    gff3_merge = makerbin + 'gff3_merge'
    cmd2 = '{1} -n -d {0}.maker.output/{0}_master_datastore_index.log'.format(seqID, gff3_merge)
    cmd2L = shlex.split(cmd2)
    p2 = subprocess.Popen(cmd2L, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p2.communicate()
    log_err('gff3_merge', seqID, stdout, stderr)

    # Generate FASTA file
    logging.info('start fasta_merge {}'.format(seqID))
    fasta_merge = makerbin + 'fasta_merge'
    cmd3 = '{1} -d {0}.maker.output/{0}_master_datastore_index.log'.format(seqID, fasta_merge)
    cmd3L = shlex.split(cmd3)
    p3 = subprocess.Popen(cmd3L, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p3.communicate()
    log_err('fasta_merge', seqID, stdout, stderr)
    
    logging.info('complete run Maker {}'.format(seqID))
    
    #delete the datastore to save disk space
    if keepStore == False:
        logging.info('remove maker output folder from MAKER working directory {}'.format(seqID))
        cmd4 = 'rm -fr {0}.maker.output'.format(seqID)
        cmd4L = shlex.split(cmd4)
        p4 = subprocess.Popen(cmd4L, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p4.communicate()
        log_err('rm MAKER datastore and blastMPI', seqID, stdout, stderr)
    
def composite_run(lock, seqID):
    logging.info('Start working on {} folder'.format(seqID))  # Chr1_2741699_2750000
    
    print(seqID)
    
    with lock:
        seqIDCoordDict = prep_maker_gff(gffdb, seqID, prefixDir)
    prep_opt_ctl(seqID, prefixDir, ctlPath)
    run_maker(seqID, prefixDir, makerBinPath, keep_datastore)
    
    return(seqIDCoordDict)
    
def run_maker_parallel(seqNameL):
    
    print(seqNameL)
    
    pool = multiprocessing.Pool(numProcess)
    manager = multiprocessing.Manager()
    lock = manager.Lock()
    seqIDCoordDictL = []
    
    try:
        seqDictCoord = pool.map(partial(composite_run, lock), seqNameL)
        seqIDCoordDictL.extend(seqDictCoord)
        pool.close()
        pool.join()
        
    except KeyboardInterrupt:
        # Allow ^C to interrupt from any thread.
        sys.stdout.write('\033[0m')
        sys.stdout.write('User Interupt\n')
    
    return(seqIDCoordDictL)

def update_gff(gffdb, prefixDir, seqIDCoord):
    '''
    output to new gff file
    :return:
    '''
    
    logging.info('[START] update gff file for {}'.format(seqIDCoord))
    
    # modify coodinate, and output modified gff file
    seqID = seqIDCoord.keys()[0]

    newGffFile = '{0}{1}/{1}.all.gff'.format(prefixDir, seqID)
    modGffFile = '{0}{1}/{1}.all.mod.gff'.format(prefixDir, seqID)
    
    # TODO if annotation file is not there, not update the file.
    
    if os.path.isfile(modGffFile):
        pass
    else:
        FILE_IN = open(newGffFile, 'rb')
        FILE_OUT = open(modGffFile, 'wb')
        adjStart = seqIDCoord[seqID]
        
        # FILE_OUT.write('##gff-version 3\n')
        
        for line in FILE_IN.readlines():
            if line.startswith('#'):
                continue
            else:
                featL = line.rstrip().split('\t')
                if featL[1] in ['maker', 'augustus_masked', 'snap_masked', 'est2genome', 'protein2genome']:
                    modCoord = [int(featL[3]) + adjStart, int(featL[4]) + adjStart + 1]
                    featL[3:5] = map(str, modCoord)

                    # featL[0] = featL[0].split('_')[0]
                    # parse seqid; e.g.: Chr9_4075176_4076008 or Scaffold_9_100_200
                    featL[0] = parse_SeqID(featL[0])[0]
                    
                    FILE_OUT.write('{0}\n'.format('\t'.join(featL)))
                    # print('\t'.join(featL))
        FILE_IN.close()
        FILE_OUT.close()
    
    logging.info('Complete update annotation db {}'.format(seqIDCoord))
    
    # load modified gff file to db
    # newGffdb = gffutils.create_db(modGffFile, dbfn=prefixDir + seqID + '.all.mod.gff.db',
    
    # logging.info('update GFF database')
    
    # db.update(modGffFile)
    
    # c = db.conn.cursor()
    # # for feature in modGffdb.features_of_type(['gene', 'mRNA']):
    # for feature in modGffdb.all_features():
    #     # print feature
    #     db._insert(feature, c)
    
    

def merge_gff(prefixDir, seqIDCoord):
    '''
    Merge modified gff files and return a merged gff file path in prefix folder
    :param seqIDCoord:
    :return:
    '''
    # merge gff file
    logging.info('[START] merge and update gff file for {}'.format(seqIDCoord))
    
    gff_fhout = open(prefixDir + 'merged_defused.all.mod.gff', 'wb')
    trans_fhout = open(prefixDir + 'merged_defused_transcripts.fa', 'wb')
    prot_fhout = open(prefixDir + 'merged_defused_protein.fa', 'wb')
    
    for seqDict in seqIDCoord:
        modGffFile = '{0}{1}/{1}.all.mod.gff'.format(prefixDir, seqDict.keys()[0])
        trans_fasta = '{0}{1}/{1}.all.maker.transcripts.fasta'.format(prefixDir, seqDict.keys()[0])
        prot_fasta = '{0}{1}/{1}.all.maker.proteins.fasta'.format(prefixDir, seqDict.keys()[0])
        gff_fn = open(modGffFile, 'rb')

        gff_fhout.write(gff_fn.read())
        gff_fn.close()
        
        if os.path.isfile(trans_fasta) and os.path.isfile(prot_fasta):
            trans_fn = open(trans_fasta, 'rb')
            prot_fn = open(prot_fasta, 'rb')
            trans_fhout.write(trans_fn.read())
            prot_fhout.write(prot_fn.read())
        else:
            continue
        
        trans_fn.close()
        prot_fn.close()

    gff_fhout.close()
    trans_fhout.close()
    prot_fhout.close()
    
    # TODO another way to do this is to merge 10 folders into one gff file, and then update the gff.db file.
    # also need to check gff has to be valid but not empty
    
    # connect to gffdb
    # try:
    #     db = gffutils.FeatureDB(gffdb)
    # except ValueError:
    #     logging.warning('database import error {}'.format(gffdb))
    #     raise ValueError('import gff/gtf sqlite data base is not present in given path')
    
    # # update gff file
    # db.update(prefixDir + 'merged_gff.all.mod.gff')
    
    logging.info('[FINISH] gff merge')
    
def write_gff(gffdb):
    
    try:
        logging.info('[START] connect gffdb to prepare output updated GFF file')
        db = gffutils.FeatureDB(gffdb)
    except ValueError:
        print(gffdb)
        raise ValueError('import gff/gtf sqlite data base is not present in given path')
    
    with open(prefixDir+'gene_features_wo_fused.gff', 'w') as fout:
        
        fout.write('## gff-version 3\n')
        
        for feature in db.all_features():
            fout.write(str(feature) + '\n')
        
        # with open(prefixDir + 'merged_gff.all.mod.gff', 'rb') as fhin:
        #     fout.write(fhin.read())
    
        # run BEDtools to sort gff file
        # logging.info('start bedtools sortbed')
        # makerExe = makerbin + 'maker'
        # cmd1 = '{maker} -q -fix_nucleotide'.format(maker=makerExe)
        # cmd1L = shlex.split(cmd1)
        # p1 = subprocess.Popen(cmd1L, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # stdout, stderr = p1.communicate()
        # log_err('MAKER', seqID, stdout, stderr)
    
    logging.info('[FINISH] write to final defused gff')
    
    
def main():
    # main routine
    flat_multiple_breaks(coordIn, 'flat_break.txt')
    seqNamesL = run_extract_seq_parallel(prefixDir + 'flat_break.txt')
    # seqNamesL = map(os.path.basename, glob.glob(prefixDir + 'Ro*')) # debugging purpose

    seqCoordL = run_maker_parallel(seqNamesL)
    del_gene_records(gffdb, coordIn)

################################################################################################
    # This part was used solely to test the gffdb update part
    #
    # seqFolderL = glob.glob(prefixDir + '/Ro*')
    # seqFolderL = map(os.path.basename, seqFolderL)
    # def parse_folder_name(folder_name):
    #     folder_list = folder_name.split('_')
    #     print(folder_list)
    #     start = folder_list[1]
    #     seq_adj = int(start) - 1
    #     return({folder_name:seq_adj})
    #
    # seqCoordL = map(parse_folder_name, seqFolderL)
################################################################################################
    
    for seqDict in seqCoordL:
        update_gff(gffdb, prefixDir, seqDict)
        
    merge_gff(prefixDir, seqCoordL)
    write_gff(gffdb)
    
    
if __name__ == '__main__':
    main()
