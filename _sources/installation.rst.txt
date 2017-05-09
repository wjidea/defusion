Installation
============
How to install `deFusion`

`deFusion` is tested with python 2.7.9.

Dependencies
------------
Install dependencies:

 #. `python/2.7.9 <http://www.python.org/>`_
 #. `gffutils/0.8.7.1 <https://github.com/daler/gffutils>`_
 #. `Biopython/1.65 <http://biopython.org/>`_
 #. `BioPerl/1.6.924 <http://bioperl.org/>`_
 #. `MAKER/r1128 <http://www.yandell-lab.org/software/maker.html>`_
 #. `BLAST+/2.2.30 <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/>`_

Prerequisites:
--------------

Under current build, `deFusion` was designed to identify and 'defuse' fused
annotated genes in MAKER annotated genomes. Output files generated from a MAKER
run are expected to be used as input files for `deFusion`:

#. Transcripts in [fasta] format
#. MAKER genome annotation in [GFF] format
#. Genome sequence in [fasta] format
#. MAKER control(ctl) files

Prepare MAKER control files
---------------------------

MAKER control files `maker_bopts.ctl`, `maker_exe.ctl`, `maker_opts.ctl` provide
the basic resource location/path informations for MAKER pipeline. Here we need
to create the MAKER controls for the deFusion pipeline, and we only need to
change the `maker_opts.ctl` file.

You can copy those three control files from your previous MAKER run, but with the
following changes:

.. role:: red

**Note: options are set to 0 unless noted in modifications**

* Section: Genome section::

    1. delete genome path

* Section: Re-annotation Using MAKER Derived GFF3::

    2. delete maker_gff path
    3. set rm_pass = 1

* Section: EST Evidence and Protein Homology Evidence::

    4. est=/path/to/your/EST/evidences.fa
    5. protein=/path/to/your/protein/evidences.fa

* Repeat Masking::

    6. set all repeat masking to null or '0'
    7. set softmask=1

* Section: Gene Prediction::

    7. snaphmm=/path/to/SNAP/snap.hmm
    8. augustus_species=augustus_spp
    9. or other gene prediction files

* Section: MAKER Behavior Options::

    10. AED_threshold=1
    11. always_complete=1
    12. keep_preds=1

Example MAKER control files: `control files <files/maker_ctl.tgz>`_


TODO
----

Using `conda <http://conda.pydata.org/docs/index.html>`_::

  conda install --channel bioconda defusion
