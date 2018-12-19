#!/bin/bash

for sra_file in *.sra; do
	#gets full path for sra file
	sra_path=$(realpath $sra_file)

	#runs fastq-dump
	fastq-dump --outdir fastq_original -Q 33 -B --split-files $sra_path
done