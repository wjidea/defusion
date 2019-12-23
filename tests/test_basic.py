from Bio.Blast.Applications import NcbiblastnCommandline

# test blastn in Bio.Blast.Applications
def blastn():
    blastfmt = "6 qseqid qlen length sstart send qstart qend evalue"
    blastnCmd = NcbiblastnCommandline(cmd='blastn', query='seq1.fa', subject='seq2.fa',
                                      evalue='1e-5', outfmt=blastfmt, max_hsps='10')
    print(blastnCmd)
    results = blastnCmd()
    print(results)

blastn()