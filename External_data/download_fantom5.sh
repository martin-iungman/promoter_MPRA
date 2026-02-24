#!/bin/bash
wget -o External_data/FANTOM5/hg38v9_samples.xlsx https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v9/basic/HumanSamples2.0.sdrf.xlsx

wget -P FANTOM5/hg38_primary_cell -r --no-parent --no-directories -A ".bed.gz" -R "index.html*" https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v9/basic/human.primary_cell.hCAGE/

wget -P FANTOM5/hg38_tissue -r --no-parent --no-directories -A ".bed.gz" -R "index.html*" https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v9/basic/human.tissue.hCAGE/

wget -P FANTOM5/hg38_cell_line -r --no-parent --no-directories -A ".bed.gz" -R "index.html*" https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v9/basic/human.cell_line.hCAGE/


