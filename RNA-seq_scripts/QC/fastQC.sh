#!/bin/bash
# This script takes all fastq files in the current directory and returns quality control figures for each, as well as a combined quality control report

## USAGE: ./fastQC.sh <directory_to_make>
# $1 <directory_to_make.txt> is the directory to be made that will contain all of the QC information and reports

## OUTPUTS:
# A variety of QC reports, but look at 'multiqc_report.html' to see a single report for all samples

## DEPENDENCIES:
# fastQC (v0.11.5)
# multiQC (v1.1)

# runs fastQC
mkdir $1
fastqc -o $1 *.fastq
cd $1
mkdir zip_files
mkdir html_files
mv *.zip zip_files
mv *.html html_files

# runs multiQC
cd zip_files
multiqc .
mv multiqc_report.html ..
cd ..
cd ..