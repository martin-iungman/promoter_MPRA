#!/bin/bash

for filename in ../FASTQ/*.fastq.gz
do 
	mv "./$filename" "./$(echo "$filename" | sed -e 's/Library_._//g')"
done
