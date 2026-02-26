#!/bin/bash
echo "#####################################################################"
echo $1
fastq1="FASTQ/Gate-$1_R1.fastq.gz"
fastq2="FASTQ/Gate-$1_R2.fastq.gz"
temp1="Trimming/TEMP_trimmed_fastq_$1_R1.fastq"
temp2="Trimming/TEMP_trimmed_fastq_$1_R2.fastq"
trimmed1="Trimming/trimmed_fastq_$1_R1.fastq"
trimmed2="Trimming/trimmed_fastq_$1_R2.fastq"
cutadapt_summary="Summary/cutadapt_fastq_$1.txt"
fastp_summary="Summary/fastp_fastq_$1.txt"
html_fastp="Summary/fastp_report_$1.html"
cutadapt -g  file:pcr_primers.fa -G file:pcr_primers.fa --discard-untrimmed -o $temp1 -p $temp2 $fastq1 $fastq2 > $cutadapt_summary
echo "Cutadapt trimming done"
fastp -i $temp1 -I $temp2 -o $trimmed1 -O $trimmed2 -q 20 --disable_adapter_trimming --cut_right --cut_right_window_size 3 --cut_right_mean_quality 20 --trim_poly_g poly_g_min_len 3 --dont_eval_duplication --html $html_fastp > $fastp_summary
echo "FASTP trimming done"
rm $temp1
rm $temp2
gzip $trimmed1
gzip $trimmed2
echo "compression done"
