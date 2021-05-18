rule clean_minimap2_cutesv:
    input:
        f"{RESULTDIR}/minimap2/cutesv/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/minimap2/cutesv/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/minimap2/cutesv/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/minimap2/cutesv/total/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/cutesv_minimap2_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmp=f"{RESULTDIR}/minimap2/cutesv/total/tmp.txt"
    shell:
        """
        echo SI00001 CUTESV > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10'| bcftools view -e 'GT[*]="RR"' | bcftools sort > {output.vcf} && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' {output.vcf} | grep -v "BND" | awk '{{FS=OFS="\t"}}{{if ($4 != "INS") print $1,$2,$3,$4; else print $1,$2,$3+$5,$4}}' | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.minimap2.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "BND" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="N") print $1,$2,$4,"BND"; else print $1,$2,$3,"BND"}}' | awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {params.tmp} && [ ! -s {params.tmp} ] && touch {output.bedpe} || pairToBed -a {params.tmp} -b <(cat {params.complex} {VCFDIR}/SI00001.minimap2.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -type neither > {output.bedpe} && rm {params.tmp} 2> {log}
        """

rule clean_ngmlr_cutesv:
    input:
        f"{RESULTDIR}/ngmlr/cutesv/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/ngmlr/cutesv/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/ngmlr/cutesv/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/ngmlr/cutesv/total/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/cutesv_ngmlr_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmp=f"{RESULTDIR}/ngmlr/cutesv/total/tmp.txt"
    shell:
        """
        echo SI00001 CUTESV > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10'| bcftools view -e 'GT[*]="RR"' | bcftools sort > {output.vcf} && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' {output.vcf} | grep -v "BND" | awk '{{FS=OFS="\t"}}{{if ($4 != "INS") print $1,$2,$3,$4; else print $1,$2,$3+$5,$4}}' | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.ngmlr.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "BND" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="N") print $1,$2,$4,"BND"; else print $1,$2,$3,"BND"}}' | awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {params.tmp} && [ ! -s {params.tmp} ] && touch {output.bedpe} || pairToBed -a {params.tmp} -b <(cat {params.complex} {VCFDIR}/SI00001.ngmlr.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -type neither > {output.bedpe} && rm {params.tmp} 2> {log}
        """

rule clean_lra_cutesv:
    input:
        f"{RESULTDIR}/lra/cutesv/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/lra/cutesv/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/lra/cutesv/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/lra/cutesv/total/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/cutesv_lra_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmp=f"{RESULTDIR}/lra/cutesv/total/tmp.txt"
    shell:
        """
        echo SI00001 CUTESV > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10'| bcftools view -e 'GT[*]="RR"' | bcftools sort > {output.vcf} && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' {output.vcf} | grep -v "BND" | awk '{{FS=OFS="\t"}}{{if ($4 != "INS") print $1,$2,$3,$4; else print $1,$2,$3+$5,$4}}' | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.lra.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "BND" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="N") print $1,$2,$4,"BND"; else print $1,$2,$3,"BND"}}' | awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {params.tmp} && [ ! -s {params.tmp} ] && touch {output.bedpe} || pairToBed -a {params.tmp} -b <(cat {params.complex} {VCFDIR}/SI00001.lra.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -type neither > {output.bedpe} && rm {params.tmp} 2> {log}
        """

rule cutesv_svs_plot:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/cutesv/total/SI00001.nobnd.bed", aligner=['minimap2', 'ngmlr', 'lra']),
        expand(f"{RESULTDIR}/{{aligner}}/cutesv/total/SI00001.bnd.bedpe", aligner=['minimap2', 'ngmlr', 'lra'])
    output:
        f"{RESULTDIR}/SI00001.cutesv.svs.pdf"
    log:
        f"{LOGDIR}/results/cutesv_svs_plot.log"
    conda: "../envs/plot.yaml"
    threads: 1
    params:
        minimap2dir = f"{RESULTDIR}/minimap2/cutesv/total",
        ngmlrdir = f"{RESULTDIR}/ngmlr/cutesv/total",
        lradir= f"{RESULTDIR}/lra/cutesv/total"
    shell:
        "Rscript {SCRIPTDIR}/circosplot.R {params.minimap2dir} {params.ngmlrdir} {params.lradir} cutesv {RESULTDIR} 2> {log}"


rule clean_minimap2_sniffles:
    input:
        f"{RESULTDIR}/minimap2/sniffles/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/minimap2/sniffles/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/minimap2/sniffles/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/minimap2/sniffles/total/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/sniffles_minimap2_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmp=f"{RESULTDIR}/minimap2/sniffles/total/tmp.txt"
    shell:
        """
        echo alignments/SI00001.minimap2.srt.bam SNIFFLES > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10' |  bcftools view -e 'GT[*]="RR"' | bcftools sort > {output.vcf} && rm {params.tmp}  && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' {output.vcf} | awk '{{FS=OFS="\t"}}{{if ($4 == "DEL" || $4 == "DUP" || $4 == "INV") print $1,$2,$3,$4; if ($4 == "INS") print $1,$2,$3+$5,$4}}' | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.minimap2.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "BND" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="N") print $1,$2,$4,"BND"; else print $1,$2,$3,"BND"}}' | awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {params.tmp} && [ ! -s {params.tmp} ] && touch {output.bedpe} ||  pairToBed -a {params.tmp} -b <(cat {params.complex} {VCFDIR}/SI00001.minimap2.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -type neither > {output.bedpe} && rm {params.tmp} 2> {log}
        """

rule clean_ngmlr_sniffles:
    input:
        f"{RESULTDIR}/ngmlr/sniffles/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/ngmlr/sniffles/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/ngmlr/sniffles/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/ngmlr/sniffles/total/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/sniffles_ngmlr_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmp=f"{RESULTDIR}/ngmlr/sniffles/total/tmp.txt"
    shell:
        """
        echo alignments/SI00001.ngmlr.srt.bam SNIFFLES > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output.vcf} && rm {params.tmp}  && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' {output.vcf} | awk '{{FS=OFS="\t"}}{{if ($4 == "DEL" || $4 == "DUP" || $4 == "INV") print $1,$2,$3,$4; if ($4 == "INS") print $1,$2,$3+$5,$4}}' | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.ngmlr.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "BND" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="N") print $1,$2,$4,"BND"; else print $1,$2,$3,"BND"}}' | awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {params.tmp} && [ ! -s {params.tmp} ] && touch {output.bedpe} ||  pairToBed -a {params.tmp} -b <(cat {params.complex} {VCFDIR}/SI00001.ngmlr.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -type neither > {output.bedpe} && rm {params.tmp} 2> {log}
        """

rule clean_lra_sniffles:
    input:
        f"{RESULTDIR}/lra/sniffles/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/lra/sniffles/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/lra/sniffles/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/lra/sniffles/total/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/sniffles_lra_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmp=f"{RESULTDIR}/lra/sniffles/total/tmp.txt"
    shell:
        """
        echo alignments/SI00001.lra.srt.bam SNIFFLES > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output.vcf} && rm {params.tmp}  && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' {output.vcf} | awk '{{FS=OFS="\t"}}{{if ($4 == "DEL" || $4 == "DUP" || $4 == "INV") print $1,$2,$3,$4; if ($4 == "INS") print $1,$2,$3+$5,$4}}' | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.lra.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "BND" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="N") print $1,$2,$4,"BND"; else print $1,$2,$3,"BND"}}' | awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {params.tmp} && [ ! -s {params.tmp} ] && touch {output.bedpe} ||  pairToBed -a {params.tmp} -b <(cat {params.complex} {VCFDIR}/SI00001.lra.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -type neither > {output.bedpe} && rm {params.tmp} 2> {log}
        """

rule sniffles_svs_plot:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/sniffles/total/SI00001.nobnd.bed", aligner=['minimap2', 'ngmlr', 'lra']),
        expand(f"{RESULTDIR}/{{aligner}}/sniffles/total/SI00001.bnd.bedpe", aligner=['minimap2', 'ngmlr', 'lra'])
    output:
        f"{RESULTDIR}/SI00001.sniffles.svs.pdf"
    log:
        f"{LOGDIR}/results/sniffles_svs_plot.log"
    conda: "../envs/plot.yaml"
    threads: 1
    params:
        minimap2dir = f"{RESULTDIR}/minimap2/sniffles/total",
        ngmlrdir = f"{RESULTDIR}/ngmlr/sniffles/total",
        lradir = f"{RESULTDIR}/lra/sniffles/total"
    shell:
        "Rscript {SCRIPTDIR}/circosplot.R {params.minimap2dir} {params.ngmlrdir} {params.lradir} sniffles {RESULTDIR} 2> {log}"

rule clean_minimap2_svim:
    input:
        f"{RESULTDIR}/minimap2/svim/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/minimap2/svim/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/minimap2/svim/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/minimap2/svim/total/SI00001.bnd.bedpe"        
    log:
        f"{LOGDIR}/results/svim_minimap2_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmp=f"{RESULTDIR}/minimap2/svim/total/tmp.txt"
    shell:
        """
        echo SI00001 SVIM > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'SUPPORT>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output.vcf} && rm {params.tmp}  && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' {output.vcf} |  grep -v "BND" | sed 's/DUP:TANDEM/DUP/g' | sed 's/DUP:INT/DUP/g' | awk '{{FS=OFS="\t"}}{{if ($4 == "DEL" || $4 == "DUP" || $4 == "INV") print $1,$2,$3,$4; else print $1,$2,$3+$5,$4}}' | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.minimap2.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "BND" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="N") print $1,$2,$4; else print $1,$2,$3}}'| awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{OFS=FS="\t"}}{{if($1!=$3) print $1"-"$2, $3"-"$4}}' | awk '{{OFS="\t"}}!seen[$1>$2 ? $1 OFS $2 : $2 OFS $1]++' | awk -F"-" '{{OFS="\t"}}$1=$1'  | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {params.tmp} && [ ! -s {params.tmp} ] && touch {output.bedpe} || pairToBed -a {params.tmp} -b <(cat {params.complex} {VCFDIR}/SI00001.minimap2.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -type neither > {output.bedpe} && rm {params.tmp} 2> {log}
        """

rule clean_ngmlr_svim:
    input:
        f"{RESULTDIR}/ngmlr/svim/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/ngmlr/svim/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/ngmlr/svim/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/ngmlr/svim/total/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/svim_ngmlr_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmp=f"{RESULTDIR}/ngmlr/svim/total/tmp.txt"
    shell:
        """
        echo SI00001 SVIM > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'SUPPORT>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output.vcf} && rm {params.tmp}  && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' {output.vcf} |  grep -v "BND" | sed 's/DUP:TANDEM/DUP/g' | sed 's/DUP:INT/DUP/g' | awk '{{FS=OFS="\t"}}{{if ($4 == "DEL" || $4 == "DUP" || $4 == "INV") print $1,$2,$3,$4; else print $1,$2,$3+$5,$4}}' | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.ngmlr.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "BND" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="N") print $1,$2,$4; else print $1,$2,$3}}'| awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{OFS=FS="\t"}}{{if($1!=$3) print $1"-"$2, $3"-"$4}}' | awk '{{OFS="\t"}}!seen[$1>$2 ? $1 OFS $2 : $2 OFS $1]++' | awk -F"-" '{{OFS="\t"}}$1=$1'  | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {params.tmp} && [ ! -s {params.tmp} ] && touch {output.bedpe} || pairToBed -a {params.tmp} -b <(cat {params.complex} {VCFDIR}/SI00001.ngmlr.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -type neither > {output.bedpe} && rm {params.tmp} 2> {log}
        """

rule clean_lra_svim:
    input:
        f"{RESULTDIR}/lra/svim/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/lra/svim/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/lra/svim/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/lra/svim/total/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/svim_lra_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmp=f"{RESULTDIR}/lra/svim/total/tmp.txt"
    shell:
        """
        echo SI00001 SVIM > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'SUPPORT>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output.vcf} && rm {params.tmp}  && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' {output.vcf} |  grep -v "BND" | sed 's/DUP:TANDEM/DUP/g' | sed 's/DUP:INT/DUP/g' | awk '{{FS=OFS="\t"}}{{if ($4 == "DEL" || $4 == "DUP" || $4 == "INV") print $1,$2,$3,$4; else print $1,$2,$3+$5,$4}}' | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.lra.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "BND" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="N") print $1,$2,$4; else print $1,$2,$3}}'| awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{OFS=FS="\t"}}{{if($1!=$3) print $1"-"$2, $3"-"$4}}' | awk '{{OFS="\t"}}!seen[$1>$2 ? $1 OFS $2 : $2 OFS $1]++' | awk -F"-" '{{OFS="\t"}}$1=$1'  | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {params.tmp} && [ ! -s {params.tmp} ] && touch {output.bedpe} || pairToBed -a {params.tmp} -b <(cat {params.complex} {VCFDIR}/SI00001.lra.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -type neither > {output.bedpe} && rm {params.tmp} 2> {log}
        """

rule svim_svs_plot:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/svim/total/SI00001.nobnd.bed", aligner=['minimap2', 'ngmlr', 'lra']),
        expand(f"{RESULTDIR}/{{aligner}}/svim/total/SI00001.bnd.bedpe", aligner=['minimap2', 'ngmlr', 'lra'])
    output:
        f"{RESULTDIR}/SI00001.svim.svs.pdf"
    log:
        f"{LOGDIR}/results/SI00001.svim_svs_plot.log"
    conda: "../envs/plot.yaml"
    threads: 1
    params:
        minimap2dir = f"{RESULTDIR}/minimap2/svim/total",
        ngmlrdir = f"{RESULTDIR}/ngmlr/svim/total",
        lradir= f"{RESULTDIR}/lra/svim/total"
    shell:
        "Rscript {SCRIPTDIR}/circosplot.R {params.minimap2dir} {params.ngmlrdir} {params.lradir} svim {RESULTDIR} 2> {log}"

rule clean_pbmm2_pbsv:
    input:
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.bnd.bedpe"        
    log:
        f"{LOGDIR}/results/pbsv_pbmm2_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmp=f"{RESULTDIR}/pbmm2/pbsv/total/tmp.txt"
    shell:
        """
        echo SI00001 PBSV > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'AD[:1]>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output.vcf} && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' {output.vcf} | grep -v "BND" | awk '{{FS=OFS="\t"}}{{if ($4 == "DEL" || $4 == "DUP" || $4 == "INV") print $1,$2,$3,$4; else print $1,$2,$3+$5,$4}}' | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.pbmm2.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "BND" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="A"||$3=="C"||$3=="G"||$3=="T") print $1,$2,$4; else print $1,$2,$3}}'| awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{OFS=FS="\t"}}{{if($1!=$3) print $1"-"$2, $3"-"$4}}' | awk '{{OFS="\t"}}!seen[$1>$2 ? $1 OFS $2 : $2 OFS $1]++' | awk -F"-" '{{OFS="\t"}}$1=$1'  | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {params.tmp} && [ ! -s {params.tmp} ] && touch {output.bedpe} || pairToBed -a {params.tmp} -b <(cat {params.complex} {VCFDIR}/SI00001.pbmm2.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -type neither > {output.bedpe} && rm {params.tmp} 2> {log}
        """

rule pbsv_svs_plot:
    input:
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.nobnd.bed",
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.bnd.bedpe"
    output:
        f"{RESULTDIR}/SI00001.pbsv.svs.pdf"
    log:
        f"{LOGDIR}/results/pbsv_svs_plot.log"
    conda: "../envs/plot.yaml"
    threads: 1
    params:
        pbmm2dir = f"{RESULTDIR}/pbmm2/pbsv/total",
    shell:
        "Rscript {SCRIPTDIR}/circosplot.R {params.pbmm2dir} pbsv {RESULTDIR} 2> {log}"

rule clean_truth:
    input:
        f"{VCFTRUTH}"
    output:
        vcf=f"{VCFDIR}/SI00001.clean.vcf",
        bed=f"{VCFDIR}/SI00001.nobnd.bed",
        bedpe=f"{VCFDIR}/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/truth_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmp=f"{VCFDIR}/tmp.txt"
    shell:
        """
        echo SI00001 TRUTH > {params.tmp} && bcftools reheader -s {params.tmp} {input} | bcftools view -f PASS | grep -v -f {params.exclude}  | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output.vcf} && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\n' {output.vcf} | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b {params.complex} -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "BND" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="N") print $1,$2,$4; else print $1,$2,$3}}' | awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{OFS=FS="\t"}}{{if($1!=$3) print $1"-"$2, $3"-"$4}}' | awk '{{OFS="\t"}}!seen[$1>$2 ? $1 OFS $2 : $2 OFS $1]++' | awk -F"-" '{{OFS="\t"}}$1=$1'  | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {output.bedpe} 2> {log}
        """

rule truth_svs_plot:
    input:
        bed=f"{VCFDIR}/SI00001.nobnd.bed",
        bedpe=f"{VCFDIR}/SI00001.bnd.bedpe"
    output:
        f"{RESULTDIR}/SI00001.truth.svs.pdf"
    log:
        f"{LOGDIR}/results/truth_svs_plot.log"
    conda: "../envs/plot.yaml"
    threads: 1
    shell:
        "Rscript {SCRIPTDIR}/circosplot.R {VCFDIR} truth {RESULTDIR} 2> {log}"

rule clean_minimap2_npinv:
    input:
        f"{RESULTDIR}/minimap2/npinv/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/minimap2/npinv/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/minimap2/npinv/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/minimap2/npinv/total/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/npinv_minimap2_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmpfile=f"{RESULTDIR}/minimap2/npinv/total/tmp.txt",
        tmpvcf=f"{RESULTDIR}/minimap2/npinv/total/tmp.vcf",
        tmpbed=f"{RESULTDIR}/minimap2/npinv/total/tmp.bed"
    shell:
        """
        echo Sample NPINV > {params.tmpfile} && grep -v -f {params.exclude} {input} | bcftools view -f PASS -i 'DV>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.tmpvcf} && bcftools query -f '%CHROM\t%POS\t%INFO/END\t%CHROM\t%POS\t%INFO/END\t%ID\t,\t+\t+\tINV\n' {params.tmpvcf} > {params.tmpbed} && rm {params.tmpvcf} && SURVIVOR bedtovcf {params.tmpbed} INV {params.tmpvcf} && rm {params.tmpbed} && bcftools reheader -s {params.tmpfile} {params.tmpvcf} > {output.vcf} && rm {params.tmpfile} && rm {params.tmpvcf} && bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\n' {output.vcf} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.minimap2.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct > {output.bed} && touch {output.bedpe} 2> {log}
        """

rule clean_ngmlr_npinv:
    input:
        f"{RESULTDIR}/ngmlr/npinv/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/ngmlr/npinv/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/ngmlr/npinv/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/ngmlr/npinv/total/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/npinv_ngmlr_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmpfile=f"{RESULTDIR}/ngmlr/npinv/total/tmp.txt",
        tmpvcf=f"{RESULTDIR}/ngmlr/npinv/total/tmp.vcf",
        tmpbed=f"{RESULTDIR}/ngmlr/npinv/total/tmp.bed"
    shell:
        """
        echo Sample NPINV > {params.tmpfile} && grep -v -f {params.exclude} {input} | bcftools view -f PASS -i 'DV>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.tmpvcf} && bcftools query -f '%CHROM\t%POS\t%INFO/END\t%CHROM\t%POS\t%INFO/END\t%ID\t,\t+\t+\tINV\n' {params.tmpvcf} > {params.tmpbed} && rm {params.tmpvcf} && SURVIVOR bedtovcf {params.tmpbed} INV {params.tmpvcf} && rm {params.tmpbed} && bcftools reheader -s {params.tmpfile} {params.tmpvcf} > {output.vcf} && rm {params.tmpfile} && rm {params.tmpvcf} && bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\n' {output.vcf} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.ngmlr.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct > {output.bed} && touch {output.bedpe} 2> {log}
        """

rule clean_lra_npinv:
    input:
        f"{RESULTDIR}/lra/npinv/total/SI00001.vcf"
    output:
        vcf=f"{RESULTDIR}/lra/npinv/total/SI00001.clean.vcf",
        bed=f"{RESULTDIR}/lra/npinv/total/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/lra/npinv/total/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/npinv_lra_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        complex = config["excluderegions"],
        tmpfile=f"{RESULTDIR}/lra/npinv/total/tmp.txt",
        tmpvcf=f"{RESULTDIR}/lra/npinv/total/tmp.vcf",
        tmpbed=f"{RESULTDIR}/lra/npinv/total/tmp.bed"
    shell:
        """
        touch {output.bedpe} && echo Sample NPINV > {params.tmpfile} && grep -v -f {params.exclude} {input} | bcftools view -f PASS -i 'DV>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.tmpvcf} && bcftools query -f '%CHROM\t%POS\t%INFO/END\t%CHROM\t%POS\t%INFO/END\t%ID\t,\t+\t+\tINV\n' {params.tmpvcf} > {params.tmpbed} && rm {params.tmpvcf} && SURVIVOR bedtovcf {params.tmpbed} INV {params.tmpvcf} && rm {params.tmpbed} && bcftools reheader -s {params.tmpfile} {params.tmpvcf} > {output.vcf} && rm {params.tmpfile} && rm {params.tmpvcf} && bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\n' {output.vcf} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.lra.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa > {output.bed} 2> {log}
        """

rule npinv_svs_plot:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/npinv/total/SI00001.nobnd.bed", aligner=['minimap2', 'ngmlr', 'lra'])
    output:
        f"{RESULTDIR}/SI00001.npinv.svs.pdf"
    log:
        f"{LOGDIR}/results/SI00001.npinv_svs_plot.log"
    conda: "../envs/plot.yaml"
    threads: 1
    params:
        minimap2dir = f"{RESULTDIR}/minimap2/npinv/total",
        ngmlrdir = f"{RESULTDIR}/ngmlr/npinv/total",
        lradir=f"{RESULTDIR}/lra/npinv/total"
    shell:
        "Rscript {SCRIPTDIR}/circosplot.R {params.minimap2dir} {params.ngmlrdir} {params.lradir} npinv {RESULTDIR} 2> {log}"


rule survivor_merge_minimap2:
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "sniffles", "svim", "npinv"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf",
        f"{VCFDIR}/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/SI00001.minimap2.upset.tsv"
    log:
        f"{LOGDIR}/results/minimap2_survivor.log"
    threads: 1
    conda: "../envs/sumsvs.yaml"
    params:
        tmp=f"{RESULTDIR}/tmpminimap2.txt",
        mergedvcf=f"{RESULTDIR}/SI00001.minimap2.merged.vcf"
    shell:
        """
        ls {input} > {params.tmp} && \
        SURVIVOR merge {params.tmp} 1000 1 1 1 0 50 {params.mergedvcf} && rm {params.tmp}  && \
        bcftools view -h {params.mergedvcf} | tail -1 | awk '{{OFS=FS="\t"}}{{print $3,$10,$11,$12,$13,$14,$15,$5,$1}}' > {output}  && \
        bcftools query -f '%CHROM\t%SVTYPE[\t%GT:%PSV:%LN:%DR:%ST:%QV:%TY:%ID:%RAL:%AAL:%CO]\n' {params.mergedvcf}  | awk '{{FS=OFS="\t"}}{{print "var"NR,$3,$4,$5,$6,$7,$8,$2,$1}}' | awk '{{FS=OFS="\t"}}{{if ($2 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $2="0"; else $2="1"; if ($3 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $3="0"; else $3="1"; if ($4 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $4="0"; else $4="1"; if ($5 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $5="0"; else $5="1";if ($6 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $6="0"; else $6="1"; if ($7 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $7="0"; else $7="1"}}1' >> {output} 2>{log}
        """

rule survivor_merge_ngmlr:
    input:
        expand(f"{RESULTDIR}/ngmlr/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "sniffles", "svim", "npinv"]),
        f"{VCFDIR}/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/SI00001.ngmlr.upset.tsv"
    log:
        f"{LOGDIR}/results/ngmlr_survivor.log"
    threads: 1
    conda: "../envs/sumsvs.yaml"
    params:
        tmp=f"{RESULTDIR}/tmpngmlr.txt",
        mergedvcf=f"{RESULTDIR}/SI00001.ngmlr.merged.vcf"
    shell:
        """
        ls {input} > {params.tmp} && \
        SURVIVOR merge {params.tmp} 1000 1 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -h {params.mergedvcf} | tail -1 |  awk '{{OFS=FS="\t"}}{{print $3,$10,$11,$12,$13,$14,$5,$1}}' > {output}  && \
        bcftools query -f '%CHROM\t%SVTYPE[\t%GT:%PSV:%LN:%DR:%ST:%QV:%TY:%ID:%RAL:%AAL:%CO]\n' {params.mergedvcf} | awk '{{FS=OFS="\t"}}{{print "var"NR,$3,$4,$5,$6,$7,$2,$1}}' | awk '{{FS=OFS="\t"}}{{if ($2 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $2="0"; else $2="1"; if ($3 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $3="0"; else $3="1"; if ($4 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $4="0"; else $4="1"; if ($5 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $5="0"; else $5="1";if ($6 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $6="0"; else $6="1"}}1' >> {output} 2>{log}
        """


rule survivor_merge_lra:
    input:
        expand(f"{RESULTDIR}/lra/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "sniffles", "svim", "npinv"]),
        f"{VCFDIR}/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/SI00001.lra.upset.tsv"
    log:
        f"{LOGDIR}/results/lra_survivor.log"
    threads: 1
    conda: "../envs/sumsvs.yaml"
    params:
        tmp=f"{RESULTDIR}/tmplra.txt",
        mergedvcf=f"{RESULTDIR}/SI00001.lra.merged.vcf"
    shell:
        """
        ls {input} > {params.tmp} && \
        SURVIVOR merge {params.tmp} 1000 1 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -h {params.mergedvcf} | tail -1 |  awk '{{OFS=FS="\t"}}{{print $3,$10,$11,$12,$13,$14,$5,$1}}' > {output}  && \
        bcftools query -f '%CHROM\t%SVTYPE[\t%GT:%PSV:%LN:%DR:%ST:%QV:%TY:%ID:%RAL:%AAL:%CO]\n' {params.mergedvcf} | awk '{{FS=OFS="\t"}}{{print "var"NR,$3,$4,$5,$6,$7,$2,$1}}' | awk '{{FS=OFS="\t"}}{{if ($2 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $2="0"; else $2="1"; if ($3 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $3="0"; else $3="1"; if ($4 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $4="0"; else $4="1"; if ($5 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $5="0"; else $5="1";if ($6 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $6="0"; else $6="1"}}1' >> {output} 2>{log}
        """

rule upset_plot:
    input:
        expand(f"{RESULTDIR}/SI00001.{{aligner}}.upset.tsv", aligner=["minimap2", "ngmlr", "lra"])
    output:
        expand(f"{RESULTDIR}/SI00001.{{aligner}}.upset.pdf", aligner=["minimap2", "ngmlr", "lra"])
    log:
        f"{LOGDIR}/results/upset_plot.log"
    threads: 1
    conda: "../envs/plot.yaml"
    shell:
        "Rscript {SCRIPTDIR}/upsetplot.R {RESULTDIR} 2> {log}"


rule survivor_consensus_minimap2:
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "sniffles", "svim", "npinv"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/SI00001.minimap2.consensus.vcf"
    log:
        f"{LOGDIR}/results/survivor_consensus_minimap2.log"
    threads: 1
    conda: "../envs/sumsvs.yaml"
    params:
        tmp=f"{RESULTDIR}/tmpminimap2_.txt"
    shell:
        """
        ls {input} > {params.tmp} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {output} && rm {params.tmp} 2> {log}
        """

rule survivor_consensus_ngmlr:
    input:
        expand(f"{RESULTDIR}/ngmlr/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "sniffles", "svim", "npinv"])
    output:
        f"{RESULTDIR}/SI00001.ngmlr.consensus.vcf"
    log:
        f"{LOGDIR}/results/survivor_consensus_ngmlr.log"
    threads: 1
    conda: "../envs/sumsvs.yaml"
    params:
        tmp=f"{RESULTDIR}/tmpngmlr_.txt"
    shell:
        """
        ls {input} > {params.tmp} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {output} && rm {params.tmp} 2> {log}
        """

rule survivor_consensus_lra:
    input:
        expand(f"{RESULTDIR}/lra/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "sniffles", "svim", "npinv"])
    output:
        f"{RESULTDIR}/SI00001.lra.consensus.vcf"
    log:
        f"{LOGDIR}/results/survivor_consensus_lra.log"
    threads: 1
    conda: "../envs/sumsvs.yaml"
    params:
        tmp=f"{RESULTDIR}/tmplra_.txt"
    shell:
        """
        ls {input} > {params.tmp} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {output} && rm {params.tmp} 2> {log}
        """

rule survivor_consensus_all:
    input:
        expand(f"{RESULTDIR}/SI00001.{{aligner}}.consensus.vcf", aligner=["minimap2", "ngmlr", "lra"])
    output:
        vcf=f"{RESULTDIR}/SI00001.consensus.vcf",
        bed=f"{RESULTDIR}/SI00001.nobnd.bed",
        bedpe=f"{RESULTDIR}/SI00001.bnd.bedpe"
    log:
        f"{LOGDIR}/results/survivor_consensus_all.log"
    threads: 1
    conda: "../envs/sumsvs.yaml"
    params:
        tmp=f"{RESULTDIR}/tmp_.txt",
        complex = config["excluderegions"]
    shell:
        """
        ls {input} > {params.tmp} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {output.vcf} && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' {output.vcf} | grep -v "TRA" | awk '{{FS=OFS="\t"}}{{if ($4 == "DEL" || $4 == "DUP" || $4 == "INV") print $1,$2,$3,$4; else print $1,$2,$3+$5,$4}}' | sort -k4 -k1,1 -k2,2n  > {params.tmp} && cut -f 4 {params.tmp} | uniq | while read sv; do awk -v var=${{sv}} '{{FS=OFS="\t"}}{{if($4==var) print $1,$2,$3,$4}}' {params.tmp} | sortBed | intersectBed -a stdin -b <(cat {params.complex} {VCFDIR}/SI00001.minimap2.exclude.tsv {VCFDIR}/SI00001.ngmlr.exclude.tsv {VCFDIR}/SI00001.pbmm2.exclude.tsv {VCFDIR}/SI00001.lra.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -v -wa | mergeBed -c 4 -o distinct >> {output.bed}; done && rm {params.tmp} && \
        bcftools query -f '%CHROM\t%POS\t%ALT\t%SVTYPE\n' {output.vcf} | grep "TRA" | awk -F"[][]" '{{OFS="\t"}}$1=$1' | awk '{{OFS="\t"}}{{if($3=="N") print $1,$2,$4; else print $1,$2,$3}}' | awk -F":" '{{OFS="\t"}}$1=$1' | awk '{{OFS=FS="\t"}}{{if($1!=$3) print $1"-"$2, $3"-"$4}}' | awk '{{OFS="\t"}}!seen[$1>$2 ? $1 OFS $2 : $2 OFS $1]++' | awk -F"-" '{{OFS="\t"}}$1=$1'  | awk '{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$4+1,"BND", 0, "+", "+"}}' > {params.tmp} && [ ! -s {params.tmp} ] && touch {output.bedpe} || pairToBed -a {params.tmp} -b <(cat {params.complex} {VCFDIR}/SI00001.minimap2.exclude.tsv {VCFDIR}/SI00001.ngmlr.exclude.tsv {VCFDIR}/SI00001.pbmm2.exclude.tsv {VCFDIR}/SI00001.lra.exclude.tsv | cut -f1-3 | sortBed | mergeBed) -type neither > {output.bedpe} && rm {params.tmp} 2> {log}
        """


rule consensus_svs_plot:
    input:
        f"{RESULTDIR}/SI00001.nobnd.bed"
    output:
        f"{RESULTDIR}/SI00001.consensus.svs.pdf"
    log:
        f"{LOGDIR}/results/consensus_svs_plot.log"
    conda: "../envs/plot.yaml"
    threads: 1
    shell:
        "Rscript {SCRIPTDIR}/circosplot.R {RESULTDIR} consensus {RESULTDIR} 2> {log}"

rule cutesv_minimap2_stats:
    input:
        f"{RESULTDIR}/minimap2/cutesv/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/minimap2/cutesv/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/cutesv_minimap2_stats.log"
    threads: 1
    params:
        outstats=f"{RESULTDIR}/minimap2/cutesv/total/stats.tsv"
    shell:
        """
        truvari stats -o {params.outstats} {input} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "CUTESV", "minimap2"}}' >> {output} 2>{log}
        """

rule cutesv_ngmlr_stats:
    input:
        f"{RESULTDIR}/ngmlr/cutesv/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/cutesv_ngmlr_stats.log"
    params:
        outstats=f"{RESULTDIR}/ngmlr/cutesv/total/stats.tsv"
    threads: 1
    shell:
        """
        truvari stats -o {params.outstats} {input} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "CUTESV", "ngmlr"}}' >> {output} 2>{log}
        """


rule cutesv_lra_stats:
    input:
        f"{RESULTDIR}/lra/cutesv/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/lra/cutesv/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/cutesv_lra_stats.log"
    params:
        outstats=f"{RESULTDIR}/lra/cutesv/total/stats.tsv"
    threads: 1
    shell:
        """
        truvari stats -o {params.outstats} {input} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "CUTESV", "lra"}}' >> {output} 2>{log}
        """


rule sniffles_minimap2_stats:
    input:
        f"{RESULTDIR}/minimap2/sniffles/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/minimap2/sniffles/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/sniffles_minimap2_stats.log"
    params:
        outstats=f"{RESULTDIR}/minimap2/sniffles/total/stats.tsv",
        tmp=f"{RESULTDIR}/minimap2/sniffles/total/tmp.vcf"
    threads: 1
    shell:
        """
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE=="INS"  || SVTYPE=="DUP" || SVTYPE=="INV" | SVTYPE=="BND"' -o {params.tmp} {input} && \
        truvari stats -o {params.outstats} {params.tmp} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && rm {params.tmp} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "SNIFFLES", "minimap2"}}' >> {output} 2>{log}
        """

rule sniffles_ngmlr_stats:
    input:
        f"{RESULTDIR}/ngmlr/sniffles/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/ngmlr/sniffles/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/sniffles_ngmlr_stats.log"
    params:
        outstats=f"{RESULTDIR}/ngmlr/sniffles/total/stats.tsv",
        tmp=f"{RESULTDIR}/ngmlr/sniffles/total/tmp.vcf"
    threads: 1
    shell:
        """
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE=="INS"  || SVTYPE=="DUP" || SVTYPE=="INV" | SVTYPE=="BND"' -o {params.tmp} {input} && \
        truvari stats -o {params.outstats} {params.tmp} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && rm {params.tmp} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "SNIFFLES", "ngmlr"}}' >> {output} 2>{log}
        """

rule sniffles_lra_stats:
    input:
        f"{RESULTDIR}/lra/sniffles/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/lra/sniffles/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/sniffles_lra_stats.log"
    params:
        outstats=f"{RESULTDIR}/lra/sniffles/total/stats.tsv",
        tmp=f"{RESULTDIR}/lra/sniffles/total/tmp.vcf"
    threads: 1
    shell:
        """
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE=="INS"  || SVTYPE=="DUP" || SVTYPE=="INV" | SVTYPE=="BND"' -o {params.tmp} {input} && \
        truvari stats -o {params.outstats} {params.tmp} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && rm {params.tmp} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "SNIFFLES", "lra"}}' >> {output} 2>{log}
        """


rule pbsv_pbmm2_stats:
    input:
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/pbmm2/pbsv/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/pbsv_pbmm2_stats.log"
    params:
        outstats=f"{RESULTDIR}/pbmm2/pbsv/total/stats.tsv"
    threads: 1
    shell:
        """
        truvari stats -o {params.outstats} {input} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "PBSV", "minimap2"}}' >> {output} 2>{log}
        """

rule svim_minimap2_stats:
    input:
        f"{RESULTDIR}/minimap2/svim/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/minimap2/svim/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/svim_minimap2_stats.log"
    threads: 1
    params:
        outstats=f"{RESULTDIR}/minimap2/svim/total/stats.tsv",
        tmp=f"{RESULTDIR}/minimap2/svim/total/tmp.vcf"
    shell:
        """
        sed 's/DUP:TANDEM/DUP/g' {input} | grep -v "DUP:INT" > {params.tmp} && \
        truvari stats -o {params.outstats} {params.tmp} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && rm {params.tmp} &&\
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "SVIM", "minimap2"}}' >> {output} 2>{log}
        """

rule svim_ngmlr_stats:
    input:
        f"{RESULTDIR}/ngmlr/svim/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/ngmlr/svim/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/svim_ngmlr_stats.log"
    threads: 1
    params:
        outstats=f"{RESULTDIR}/ngmlr/svim/total/stats.tsv",
        tmp=f"{RESULTDIR}/ngmlr/svim/total/tmp.vcf"
    shell:
        """
        sed 's/DUP:TANDEM/DUP/g' {input} | grep -v "DUP:INT" > {params.tmp} && \
        truvari stats -o {params.outstats} {params.tmp} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && rm {params.tmp} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "SVIM", "ngmlr"}}' >> {output} 2>{log}
        """

rule svim_lra_stats:
    input:
        f"{RESULTDIR}/lra/svim/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/lra/svim/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/svim_lra_stats.log"
    threads: 1
    params:
        outstats=f"{RESULTDIR}/lra/svim/total/stats.tsv",
        tmp=f"{RESULTDIR}/lra/svim/total/tmp.vcf"
    shell:
        """
        sed 's/DUP:TANDEM/DUP/g' {input} | grep -v "DUP:INT" > {params.tmp} && \
        truvari stats -o {params.outstats} {params.tmp} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && rm {params.tmp} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "SVIM", "lra"}}' >> {output} 2>{log}
        """

rule npinv_minimap2_stats:
    input:
        f"{RESULTDIR}/minimap2/npinv/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/minimap2/npinv/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/npinv_minimap2_stats.log"
    threads: 1
    params:
        outstats=f"{RESULTDIR}/minimap2/npinv/total/stats.tsv"
    shell:
        """
        truvari stats -o {params.outstats} {input} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "NPINV", "minimap2"}}' >> {output} 2>{log}
        """

rule npinv_ngmlr_stats:
    input:
        f"{RESULTDIR}/ngmlr/npinv/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/ngmlr/npinv/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/npinv_ngmlr_stats.log"
    threads: 1
    params:
        outstats=f"{RESULTDIR}/ngmlr/npinv/total/stats.tsv"
    shell:
        """
        truvari stats -o {params.outstats} {input} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "NPINV", "ngmlr"}}' >> {output} 2>{log}
        """

rule npinv_lra_stats:
    input:
        f"{RESULTDIR}/lra/npinv/total/SI00001.clean.vcf"
    output:
        f"{RESULTDIR}/lra/npinv/total/svdist.tsv"
    log:
        f"{LOGDIR}/results/npinv_lra_stats.log"
    threads: 1
    params:
        outstats=f"{RESULTDIR}/lra/npinv/total/stats.tsv"
    shell:
        """
        truvari stats -o {params.outstats} {input} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "NPINV", "lra"}}' >> {output} 2>{log}
        """

rule truth_stats:
    input:
        f"{VCFDIR}/SI00001.clean.vcf"
    output:
        f"{VCFDIR}/svdist.tsv"
    log:
        f"{LOGDIR}/results/truth_stats.log"
    threads: 1
    params:
        outstats=f"{VCFDIR}/stats.tsv"
    shell:
        """
        truvari stats -o {params.outstats} {input} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "TRUTH", "truth"}}' >> {output} 2>{log}
        """

rule combine_stats:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/{{tool}}/total/svdist.tsv", aligner=["minimap2", "ngmlr", "lra"], tool=["cutesv", "sniffles", "svim", "npinv"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/svdist.tsv",
        f"{VCFDIR}/svdist.tsv"
    output:
        f"{RESULTDIR}/SI00001.varstats.tsv"
    log:
        f"{LOGDIR}/results/combine_stats.log"
    threads: 1
    shell:
        """
        awk "FNR>1 || NR==1" {input} > {output} 2>{log}
        """





