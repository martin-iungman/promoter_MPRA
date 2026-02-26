#!/bin/bash
bam_file="../Alignment/alignment_$1.bam"
sorted_file="../Alignment/alignment_sorted_$1.bam"
samtools sort -o $sorted_file $bam_file
samtools index $sorted_file
rm $bam_file