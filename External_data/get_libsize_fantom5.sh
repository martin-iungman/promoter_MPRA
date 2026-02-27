#!/bin/bash
# the input should be the name of the desired folder of FANTOM (tissue, cell line or primary cell)
cd FANTOM5/$1
ls | cat > ../$1_files.txt
n=(ls -l | grep -v ^d | wc -l)
touch ../$1_libsize.txt
for i in $(seq 1 $n)
do
a="${i}p"
file_gz=$(sed -n $a ../$1_files.txt)
file_bed=$(echo $file_gz | sed 's/.gz//')
gzip -d $file_gz
size=$(cut -d' ' -f5 $file_bed | cut -f5 | awk '{s+=$1} END {print s}')
echo "${size} ${file_bed}" | cat >> ../$1_libsize.txt
gzip $file_bed
done
