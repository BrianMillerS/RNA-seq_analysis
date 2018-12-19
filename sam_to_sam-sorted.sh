#!/bin/bash
# sorts all sam files in curdir
# sorted by coordinates

for sam_file in *.sam; do
	#takes the file name root from each sam_file (e.g. hello.sam -> hello_sorted.sam)
	output_file=$(echo $sam_file | cut -d'.' -f 1)
	output_file+="_sorted_pos.bam"

	#runs picard
	java -jar $PICARD SortSam I=$sam_file O=$output_file SORT_ORDER=coordinate
done


