Tutorial
========

.. image:: images/diagram.png

Prepare your prerequisites
--------------------------

1. install dependencies
2. prepare outputs from previous MAKER run
3. run detect_fused_gene.py
4. run extract_to_maker.py
5. review results or trouble shooting


Example - rice chromosome 9
---------------------------

1. detect fused gene candidates::

	python bin/1_detect_fused_gene.py -i resources/rice_chr9_transcript.fa -g resources/rice_chr9.gff -n 30 -p chr9_S1R1

2. split fused gene and locally re-annotation::

	python bin/2_extract_to_maker.py -i resources/rice_chr9.fa -d chr9_S1R1/gff.db -c chr9_S1R1/break_coordinates.brk \
	-t resources/rice_chr9_ctl/ -n 28 -p chr9_S2R1

3. drop fused gene entries from original transcript and protein sequence fasta files::
	
	mkdir chr9_post_process
	
	python bin/3_process_fasta_file.py -i resources/rice_chr9_transcript.fa -g chr9_S2R1/gene_features_wo_fused.gff  \
	-o chr9_post_process/rice_chr9_transcript_drop_fused.fa
	
	python bin/3_process_fasta_file.py -i resources/rice_chr9_protein.fa -g chr9_S2R1/gene_features_wo_fused.gff  \
	-o chr9_post_process/rice_chr9_protein_drop_fused.fa

4. run maker standard on defused gene set. MAKER standard is a procedure to get rid of low-quality gene models::
	
	python bin/4_run_maker_standard.py -t chr9_S2R1/merged_defused_transcripts.fa -p chr9_S2R1/merged_defused_protein.fa \
	-g chr9_S2R1/merged_defused.all.mod.gff -a Pfam/Pfam-A.hmm -o chr9_post_process/
	

5. generate AED score improvement plot::
	
	python bin/5_generate_report.py -b chr9_S1R1/break_coordinates.brk \
	-i resources/rice_chr9.gff \
	-g chr9_post_process/merged_defused.all.mod.std.gff \
	-o chr9_post_process/chr9


6. final assembly::
	
	mv chr9_S2R1/gene_features_wo_fused.gff chr9_post_process/
	
	cd chr9_post_process
	
	cat gene_features_wo_fused.gff merged_defused.all.mod.std.gff > chr9_V1.gff
	cat rice_chr9_transcript_drop_fused.fa merged_defused_transcripts.std.fa > chr9_V1_transcripts.fa
	cat rice_chr9_protein_drop_fused.fa merged_defused_protein.std.fa > chr9_V1_protein.fa
	
7. DONE Final outputs::
	
	chr9_V1.gff
	chr9_V1_transcripts.fa
	chr9_V1_protein.fa
	
