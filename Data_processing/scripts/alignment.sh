#!/bin/bash

trimmed1="../Trimming/trimmed_fastq_$1_R1.fastq.gz"
trimmed2="../Trimming/trimmed_fastq_$1_R2.fastq.gz"
sam_file="../Alignment/alignment_$1.sam"
hisat_summary_file="../Summary/alignment_summary_$1.txt"
hisat_summary_file2="../Summary/hisat2_summary_$1.txt"
bam_file="../Alignment/alignment_$1.bam"
hisat2 --no-spliced-alignment -I 230 -X 270 --no-discordant -x ../library_index/library_idx -1 $trimmed1 -2 $trimmed2 -S $sam_file --summary-file $hisat_summary_file
samtools view -b -F 256 $sam_file > $bam_file
gzip $sam_file
