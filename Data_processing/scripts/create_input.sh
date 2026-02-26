#!/bin/bash
for Gate in 1 2 3 4 5 6
do
	for Lib in 1 2
	do
		echo ${Gate}-${Lib} >> tmp.txt
	done
done
