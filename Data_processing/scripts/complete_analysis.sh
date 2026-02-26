#!/bin/bash

source create_input.sh
#source rename_fastq.sh

mkdir ../FASTQC
mkdir ../Trimming
mkdir ../Alignment
mkdir ../Count_table
mkdir ../Summary

echo "Start QC"
fastqc -t12 --outdir ../FASTQC ../FASTQ/*.fastq
echo "QC done"

parallel -j12 source trimming.sh :::: tmp.txt
echo "Trimming done"

echo "Indexing library"
hisat2-build ../../Library_data/data/promoter_wo_dupl_SpikeIns.fa library_index/library_idx

parallel -j12 source alignment.sh :::: tmp.txt
echo "alignment done"

parallel -j12 source bam_process.sh :::: tmp.txt
echo "BAM files sorted and indexed"

parallel -j12 source count_reads.sh :::: tmp.txt
echo "Count matrix made"

rm  tmp.txt