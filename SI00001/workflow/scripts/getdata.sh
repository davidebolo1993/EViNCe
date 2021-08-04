#!/bin/bash

cd ../../data

# get GRCh38 reference with decoys

while true; do

	wget -T 15 -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa && break

done

samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

#get chromosomes that we need to exclude later

cut -f1 GRCh38_full_analysis_set_plus_decoy_hla.fa.fai | tail -n +25 > excludechroms.txt

#simulate

VISOR HACk -g GRCh38_full_analysis_set_plus_decoy_hla.fa -b hack.bed -o hack_out && cp GRCh38_full_analysis_set_plus_decoy_hla.fa hack_out/h2.fa && samtools faidx hack_out/h2.fa
cut -f1,2 <(head -24 hack_out/h1.fa.fai) <(head -24 hack_out/h2.fa.fai) <(head -24 GRCh38_full_analysis_set_plus_decoy_hla.fa.fai) > haplochroms.dim.tsv
sort haplochroms.dim.tsv | awk '$2 > maxvals[$1] {lines[$1]=$0; maxvals[$1]=$2} END { for (tag in lines) print lines[tag] }' | sort -k 1,1 -V -s > maxdims.tsv
awk 'OFS=FS="\t"''{print $1, "1", $2, "100.0", "100.0"}' maxdims.tsv > laser.bed && head -22 laser.bed > tmp.bed && tail -n+23 laser.bed | sed 's,\t100.0\t,\t50.0\t,g' >> tmp.bed && mv tmp.bed laser.bed
rm haplochroms.dim.tsv maxdims.tsv
VISOR LASeR -g GRCh38_full_analysis_set_plus_decoy_hla.fa -b laser.bed -s hack_out -o laser_out --threads 30 --coverage 50 --identity_min 88 --identity_max 98 --identity_stdev 5 --junk_reads 1 --random_reads 1 --chimera_reads 1 --glitches_rate 10000 --glitches_size 25 --glitches_skip 25
samtools fastq -@ 30 laser_out/sim.srt.bam > SI00001.fastq && rm -r laser_out && find "$(pwd)" -name "*.fastq" > input.fofn

cd - 
