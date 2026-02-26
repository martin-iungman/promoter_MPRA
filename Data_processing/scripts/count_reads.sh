#!/bin/bash

bam_file="Alignment/alignment_sorted_$1.bam"
output="Count_table/counts_ampliconseq2023_$1.tsv"	
for seq in $(cat ../Library_data/res/seq_id.txt); do
	num_reads=$(samtools view -c -f 2 $bam_file $seq)
	echo -e "$seq\t$num_reads"
	done > $output
