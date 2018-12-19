#!/bin/bash
#This script converts all .sam files in current directory to .bam files (sorted by coordinates) and makes a .bai file for each bam file 

## DEPENDENCIES:
#PICARD (2.9.4)
#samtools (0.1.19)

#PICARD="$HOME/opt/picard-2.9.4.jar"

# make bam files
for sam_file in *.sam; do
	#takes the file name root from each sam_file (e.g. hello.sam -> hello.bam)
	bam_file=$(echo $sam_file | cut -d'.' -f 1)
	bam_file+=".bam"

	#runs picard
	java -jar $PICARD SortSam I=$sam_file O=$bam_file SORT_ORDER=coordinate
done

# make bai (bam indexed files)
for bam_file in *.bam; do
	samtools indsex $bam_file
done
