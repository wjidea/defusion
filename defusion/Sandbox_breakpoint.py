import sys
import logging
import time
import gffutils

startTime = time.time()

def breakPont():
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
    init = 0
    # determine the breaking point TODO improve breaking point selection composite score
    # with open(blastOut)
    
    # implement the simple version of this function divide them in the mid point
    dbName = '/Volumes/JW_REGULAR/tmp_rice/gff.db'
    
    try:
        db = gffutils.FeatureDB(dbName)
        logging.debug('Found gff.db, now loading')
        # db = sqlite3.connect(dbName)
    except ValueError:
        print('Check gff.db is not in the working directory')
        sys.exit()
    
    # find the mid-point and fit into the mid of the intron
    def inRange(point, interval):
        if point >= min(interval) and point <= max(interval):
            return (True)
        else:
            return (False)
    
    # maker-Chr10-snap-gene-129.36-mRNA-1   Chr10:12968107..12983848
    # maker-Chr4-snap-gene-218.40-mRNA-1      Chr4:21842221..21864472 15
    # maker-Chr7-augustus-gene-258.26-mRNA-1  Chr7:25886858..25887977 2
    # augustus_masked-Chr2-processed-gene-29.6-mRNA-1 Chr2:2945613..2955585   3
    # maker-Chr1-snap-gene-27.45-mRNA-1       Chr1:2741699..2761866   17
    
    # db.execute('''
    #     CREATE INDEX chr_end_start_index ON Features(seqid, start, end)
    # ''')
    
    ROI = db.region(seqid='Chr10', start=12968107, end=12983848, completely_within=False)
    gene = db['maker-Chr1-snap-gene-27.45-mRNA-1']
    

    featKept = []
    
    
    transcriptLen = 0
    for i in ROI:
        if i.featuretype == 'gene' and i.strand == gene.strand:
            print(i.source, i.featuretype, i.strand, len(i), i.start, i.end)
            featKept.append(i)
        if i.featuretype == 'exon' and i.strand == gene.strand:
            transcriptLen += len(i)
            print(i.source, i.featuretype, i.strand, len(i), i.start, i.end, transcriptLen)
        if i.source == 'augustus_masked' and i.featuretype == 'match' and i.strand == gene.strand:
            print(i.source, i.featuretype, i.strand, len(i), i.start, i.end)
        if i.source == 'snap_masked' and i.featuretype == 'match' and i.strand == gene.strand:
            print(i.source, i.featuretype, i.strand, len(i), i.start, i.end)
    
    for j in featKept:
        print(j)

breakPont()

timePassed = time.time() - startTime

print(timePassed)