# RNA-seq_scripts
Scripts for mRNA-seq: fetching SRR data, quality control, alignment, counting, differential expression, and gene ontology analysis

These sets of scripts are to be used at the various stages of RNA-Seq analysis for detection of differential expression followed by a gene ontology analysis.

File Descriptions:
step1_SRR_to_fq.sh (convert downloaded SRR files to fastq files)
step2_fastQC.sh (generate a quality control report for multiple fastq files)
step2_trimmomatic.sh (perform quality control steps on multiple fastq files)
step3_alignment_and_counting.sh (align fastq reads to a genome and count the alignments)
step4_differerntial_expression_and_gene_ontology_analysis.R (compare counts results and calculate significantly differentially expressed genes and significantly enriched gene ontologies)

Other scripts (used to convert SAM files to BAM files):
sam_to_bam.sh
sam_to_sam-sorted.sh