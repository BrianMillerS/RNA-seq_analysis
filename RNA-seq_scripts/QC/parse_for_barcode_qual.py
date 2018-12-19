#!/usr/bin/python2.7
import time
from Bio import SeqIO
import sys


"""
This script reads in a single fastq file and returns a single fastq file excluding reads that have any of the barcode nucleotides
(first three) below a quality threshold.
"""
start_time = time.time()
print("\nEntering python script: 'parse_for_barcode_qual.py'")

# inputs
fastq_file = sys.argv[1]  # fastq file to parse
qual_threshold = int(sys.argv[2])  # int for the quality desired
output_file = sys.argv[3]  # file to be written

filtered_fastq = open(sys.argv[3], "a")  # open results file
# read your fastq file
for read in SeqIO.parse(fastq_file, 'fastq'):
    if read.letter_annotations["phred_quality"][0] < qual_threshold:  # check quality at position 0
        continue
    elif read.letter_annotations["phred_quality"][1] < qual_threshold:  # check quality at position 1
        continue
    elif read.letter_annotations["phred_quality"][2] < qual_threshold:  # check quality at position 2
        continue
    else:
        SeqIO.write(read, filtered_fastq, 'fastq')  # write the record if all qualities are above qual_threshold

# Close file
filtered_fastq.close()

print("...time elapsed: " + str(time.time() - start_time) + " seconds")
print("...exiting python script: 'parse_for_barcode_qual.py'")
