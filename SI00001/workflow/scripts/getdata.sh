#!/bin/bash

cd ../../data

# get GRCh38 reference with decoys

while true; do

	wget -T 15 -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa && break

done

samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

#get chromosomes that we need to exclude later

cut -f1 GRCh38_full_analysis_set_plus_decoy_hla.fa.fai | tail -n +25 > excludechroms.txt

#get regions we do not want to include in final calls (telomeres/centromeres/lowcomplexity)

wget https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg38.excl.tsv && head -172 human.hg38.excl.tsv | grep "^chr"  > excluderegions.txt
sort -k1,1 -k2,2n excluderegions.txt > excluderegions.srt.bed

rm human.hg38.excl.tsv excluderegions.txt


cd -
