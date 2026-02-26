#!/bin/bash

input1="Trimming/Trimmed/trimmed_fastq_$1_R1.fastq"
input2="Trimming/Trimmed/trimmed_fastq_$1_R2.fastq"
output1="FASTP/fastp_trimmed_Gate-$1_R1.fastq"
output2="FASTP/fastp_trimmed_Gate-$1_R2.fastq"
unpaired1="FASTP/fastp_unpaired_trimmed_Gate-$1_R1.fastq"
unpaired2="FASTP/fastp_unpaired_trimmed_Gate-$1_R2.fastq"
html="FASTP/fastp_$1.html"
fastp -i $input1 -I $input2 -o $output1 -O $output2 --unpaired1 $unpaired1 --unpaired2 $unpaired2 -q 20 --disable_adapter_trimming --cut_right --cut_right_window_size 3 --cut_right_mean_quality 20 --trim_poly_g poly_g_min_len 3 --dont_eval_duplication --html $html
