#!/bin/bash
# the input should be the name of the desired folder of FANTOM (tissue, cell line or primary cell)

cd FANTOM5
mkdir Library_filt_CAGE/
mkdir Library_filt_CAGE/${1}
seqs=$(wc -l < ${1}_files.txt

for i in $(seq 1 $seqs)
do
input=$(sed -n ${i}p {$1}_files.txt)
bedtools intersect -a "$1/$input" -b ../../Library_data/res/library.bed > "Library_filt_CAGE/$1/${input}" 
done
