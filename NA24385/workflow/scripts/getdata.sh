#!/bin/bash

#get NA24385 ultra-long PromethION calls on hs37d5

mkdir ../../data 
cd ../../data

#get calls for NA24385

while true;do

	wget -T 15 -c  ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed && break

done

while true; do

	wget -T 15 -c ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz && break

done

#get NA24385 ultra-long PromethION raw data

while true; do
	
	wget -T 15 -c ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_1.fastq.gz && break

done

while true; do

	wget -T 15 -c ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_2.fastq.gz && break

done

while true; do

	wget -T 15 -c ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_3.fastq.gz && break

done

# get hs37d5 reference (hg19 with decoys)

while true; do

	wget -T 15 -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz && break

done

#also decompress reference as some tools does not support .gz format

gunzip -c hs37d5.fa.gz > hs37d5.fa

#also create .fofn with fastq input for pbmmi aligner

find "$(pwd)" -name "*.fastq.gz" > input.fofn

#index VCF file

tabix HG002_SVs_Tier1_v0.6.vcf.gz

#get chromosomes that we need to exclude later

cut -f1 hs37d5.fa.fai | tail -n +25 > excludechroms.txt

#get regions we do not want to include in final calls (telomeres/centromeres/lowcomplexity/pericentromeres/subtelomeres)

wget https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg19.excl.tsv && head -144 human.hg19.excl.tsv | grep -v "^chr"  > excluderegions.txt
#wget https://raw.githubusercontent.com/shenlab-sinai/diffreps/master/lib/PJ/Database/hg19/hg19.pericentromeres.bed && awk 'OFS=FS="\t"''{print $0, "pericentromere"}' hg19.pericentromeres.bed | sed  's/^chr//g' >> excluderegions.txt
#wget https://raw.githubusercontent.com/shenlab-sinai/diffreps/master/lib/PJ/Database/hg19/hg19.subtelomeres.bed && awk 'OFS=FS="\t"''{print $0, "subtelomeres"}' hg19.subtelomeres.bed | sed  's/^chr//g' >> excluderegions.txt
sort -k1,1 -k2,2n excluderegions.txt > excluderegions.srt.bed
rm human.hg19.excl.tsv excluderegions.txt #hg19.pericentromeres.bed hg19.subtelomeres.bed


cd -

