#!/bin/sh
#Genome (just canonical chromosomes)
wget -O External_data/genome/GCF_000001405.40_GRCh38.p14_genomic.fna https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED 
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' External_data/genome/GCF_000001405.40_GRCh38.p14_genomic.fna | grep -E "^>NC" | tr "\t" "\n" > External_data/genome/GRCh38.p14_canonical.fa

#CpG islands
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz    | gunzip -c    | awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4, $5$6, $7, $8, $9, $10, $11, $12 }'   > External_data/cpgIslandExt.hg38.bed

##PhyloP score
wget -O External_data/hg38.phyloP100way.bw https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw 

#Functional conservation
wget -O External_data/Sup1_mouse_human_Young2015.txt https://genome.cshlp.org/content/suppl/2015/08/17/gr.190546.115.DC1/Supplemental_File1.txt

#sample ontology
wget  -O External_data/Tissue_cell_ontology https://static-content.springer.com/esm/art%3A10.1038%2Fnature12787/MediaObjects/41586_2014_BFnature12787_MOESM8_ESM.zip | gzip 

#Repeated Elements
# change with nano first three lines for: "SW_score perc_div perc_del perc_ins seqname start end genome_left cons_strand matching_repeat repeat_class repeat_start repeat_end repeat_left ID",
wget https://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz | gzip -c | sed -i '1,2,3d' | sed -r -i 's/^\s+//g' 
echo "SW_score perc_div perc_del perc_ins seqname start end genome_left cons_strand matching_repeat repeat_class repeat_start repeat_end repeat_left ID" > tmp.txt
cat tmp.txt hg38.fa.out > hg38.fa.out
rm tmp.txt

# FANTOM5 enhancer
wget -O External_data/F5.hg38.enhancers.bed https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz

#ChIP-atlas
wget https://chip-atlas.dbcls.jp/data/hg38/allPeaks_light/allPeaks_light.hg38.50.bed.gz
wget https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab
bedtools intersect -wao -a ../Library_data/res/library.bed -b allPeaks_light.hg38.50.bed.gz > allPeaks_light.hg38.50_lib.bed
# sort & cut
awk -F'\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' experimentList.tab |  awk -F'\t' '$2=="hg38"'| sort -k1 > experimentList_chipatlas.tab
sort -k10  allPeaks_light.hg38.50_lib.bed > allPeaks_light.hg38.50_lib_sorted.bed
join -1 10 -2 1 -t $'\t' allPeaks_light.hg38.50_lib_sorted.bed experimentList_chipatlas.tab > allPeaks_light.hg38.50_lib_joined.bed

#PUFFIN 
wget -O External_data/dudnyk2024_sup.zip https://www.science.org/doi/suppl/10.1126/science.adj0116/suppl_file/science.adj0116_tables_s1_to_s8.zip | gzip -d

#remap
wget https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_nr_macs2_hg38_v1_0.bed.gz | gzip -d

#EPD
wget -O External_data/human38_epdnew.bed https://epd.expasy.org/ftp/epdnew/H_sapiens/current/Hs_EPDnew.bed