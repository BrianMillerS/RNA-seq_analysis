#!/bin/bash
# This script runs hisat2 and featurecounts on a set of fastq files to generate counts files

## USAGE: ./Full_Alignment_and_Counting.sh <counts_output_reverse_filename.txt> <counts_output_forward_filename.txt> <genome.fasta> <annotation_file.gtf> <sam_add-on>
# $1 <counts_output_reverse_filename.txt> is the name of a file that will be created by featureCounts containing the counts for all samples for every gene using the -s 2 parameter
# $2 <counts_output_forward_filename.txt> is the name of a file that will be created by featureCounts containing the counts for all samples for every gene using the -s 1 parameter
# $3 <genome.fasta> is the genome filename to be used
# $4 <annotation_file.gtf> is the annotation filename (in gtf format) to be used
# $5 <sam_add-on> is the word that you want added to the end of the .sam filenames (e.g. if <sam add-on> = "Couvillion" then all sam files will have a *Couvillion.sam ending)

## INPUTS:
# - refrence genome in .fasta format
# - anotation file in .gtf format
# - multiple fastq files

## OUTPUTS:
# Alignment files in .sam format
# A featurecounts counts txt file with reverse counts for every annotaion (gene if transcriptome)
# A featurecounts counts txt file with forward counts for every annotaion (gene if transcriptome)

## DEPENDENCIES:
# hisat2 (2.1.0)
# featurecounts (v1.5.3)
#
# Test given lines below individually in Terminal, programs should be already added to your $PATH variable:
# hisat2
# hisat2_extract_splice_sites.py
# hisat2_extract_exons.py
# hisat2-build
# featureCounts

cwd=$(pwd)
cwd+="/"
echo "Directory: $cwd"
echo "Date: $(date)"
echo
echo "STEP 1: Creating pre files needed for alignment(sam)-index creation"
hisat2_extract_splice_sites.py $4 >>splice_sites.txt
hisat2_extract_exons.py $4 >>exons.txt
echo
echo "STEP 2: Creating index files"
hisat2-build -f --ss splice_sites.txt --exon exons.txt $3 ht2_ref
echo
echo "STEP 3: Running Hisat2 on each fastq file, output is a .sam file for each fastq"
for fastq_file in /media/pollardlab/POLLARDLAB3/RNAi_Project_T_thermophila_RNA-seq/Analysis_V2/STEP1_Original_files_and_QC/31_barcode_quality/trimmCleaned_fastq_files/*.fastq; do
	file_name=$(echo $fastq_file | cut -d'.' -f 1)  # takes the file name root from each fastq_file (e.g. bob/hello.fastq -> bob/hello)
	file_name=$(basename $file_name)  # takes the file name root from each ffile and excludes the rest of the path (e.g. bob/hello -> hello)
	echo
	echo "OUTPUT for $file_name"
	echo
	file_name+="_$5.sam"  # adds .sam to the filename

	#runs hisat2, aligning to the forward strand
	hisat2 -k 1 -q --phred33 --known-splicesite-infile splice_sites.txt --time --rna-strandness "R" -x ht2_ref -U $fastq_file -S $cwd$file_name 
	# explaining the above line of code:
	# [-k 1] find the single best alignment for a given read
	# [-q] tells hisat2 that the input file is a fastq
	# [--phred33] tells hisat2 that the fastq is quality phred 33 encoded
	# [--known-splicesite-infile splice_sites.txt] gives hisat2 the slice sites file
	# [--time] prints walltime taken for a given alignment
	# [--rna-strandness "R"] tells hisat2 to align to the reverse strand
	# [-x ht2_ref] gives hisat2 the name of the exon locations file
	# [-U $fastq_file] input fastq
	# [-S $file_name ] output sam file name
done


#move files into folders
mkdir alignment_refrence_files
mv exons.txt splice_sites.txt *.ht2 alignment_refrence_files
mkdir sam_files
cp $4 sam_files  # copy gtf into the sam_files folder, needed for step 4
mv *.sam sam_files

#Generate count file from all of the .sam files
echo "STEP 4: Running featureCounts using [-s 2] on each sam file, output is the single txt counts file $1"
cd sam_files
featureCounts -s 2 -t exon -g gene_name -a $4 -o $1 *.sam
# explaining the above line of code:
# [-s 2] tells featureCounts to count on the reverse strand (NOTE this is not redundant to the "R" during)
# [-t exon] tells featureCounts the feature type that is being counted 
# [-g gene_name] tells featureCounts to report the counts using this name identifier from the gtf file (in every case it is 'gene_name' eventhough our custom annotations are not genes)
# [-a $4] gives featureCounts the annotation set, gtf file
# [-o $1] tells featureCounts what the counts file name should be
# [*.sam] tells featureCounts to run on all .sam files in the current working directory

#Generate count file from all of the .sam files
echo "STEP 5: Running featureCounts using [-s 1] on each sam file, output is the single txt counts file $2"
cd sam_files
featureCounts -s 1 -t exon -g gene_name -a $4 -o $2 *.sam
# explaining the above line of code:
# [-s 1] tells featureCounts to count on the forward strand
# [-t exon] tells featureCounts the feature type that is being counted 
# [-g gene_name] tells featureCounts to report the counts using this name identifier from the gtf file (in every case it is 'gene_name' eventhough our custom annotations are not genes)
# [-a $4] gives featureCounts the annotation set, gtf file
# [-o $2] tells featureCounts what the counts file name should be
# [*.sam] tells featureCounts to run on all .sam files in the current working directory

cd ..
echo "Date: $(date)"
echo "Full_Alignment.sh done."