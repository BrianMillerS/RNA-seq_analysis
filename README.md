# RNA-seq_analysis
Scripts for mRNA-seq: fetching SRR data, quality control, alignment, counting, differential expression, and gene ontology analysis

These sets of scripts are to be used at the various stages of RNA-Seq analysis for detection of differential expression followed by a gene ontology analysis.

File Descriptions:<br/>
step1_SRR_to_fq.sh (uses fastq-dump to convert downloaded SRR files to fastq files)<br/>
step2_fastQC.sh (uses fastqc and multiqc to generate a quality control report for multiple fastq files)<br/>
step2_trimmomatic.sh (uses trimmomatic to perform quality control steps on multiple fastq files)<br/>
step3_alignment_and_counting.sh (uses hisat2 and featureCounts to align fastq reads to a genome and count the alignments)<br/>
step4_differerntial_expression_and_gene_ontology_analysis.R (uses Deseq2, GOstats, and various other packages to compare counts results and determine significantly differentially expressed genes and significantly enriched gene ontologies from sets of genes)<br/>

Other scripts (used to convert SAM files to BAM files):<br/>
sam_to_bam.sh<br/>
sam_to_sam-sorted.sh<br/>
