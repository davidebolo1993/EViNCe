rule clean_minimap2_cutesv:
    input:
        f"{RESULTDIR}/minimap2/cutesv/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/minimap2/cutesv/total/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/cutesv_minimap2_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/minimap2/cutesv/total/tmp.txt"
    shell:
        """
        echo SI00001 CUTESV > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10'| bcftools view -e 'GT[*]="RR"' | bcftools sort > {output} && rm {params.tmp} 2> {log}
        """

rule clean_ngmlr_cutesv:
    input:
        f"{RESULTDIR}/ngmlr/cutesv/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/total/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/cutesv_ngmlr_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/ngmlr/cutesv/total/tmp.txt"
    shell:
        """
        echo SI00001 CUTESV > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10'| bcftools view -e 'GT[*]="RR"' | bcftools sort > {output} && rm {params.tmp} 2> {log}
        """

rule clean_lra_cutesv:
    input:
        f"{RESULTDIR}/lra/cutesv/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/lra/cutesv/total/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/cutesv_lra_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/lra/cutesv/total/tmp.txt"
    shell:
        """
        echo SI00001 CUTESV > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10'| bcftools view -e 'GT[*]="RR"' | bcftools sort > {output} && rm {params.tmp} 2> {log}
        """

rule clean_minimap2_sniffles:
    input:
        f"{RESULTDIR}/minimap2/sniffles/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/minimap2/sniffles/total/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/sniffles_minimap2_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/minimap2/sniffles/total/tmp.txt"
    shell:
        """
        echo alignments/SI00001.minimap2.srt.bam SNIFFLES > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10' |  bcftools view -e 'GT[*]="RR"' | bcftools sort > {output} && rm {params.tmp} 2> {log}
        """

rule clean_ngmlr_sniffles:
    input:
        f"{RESULTDIR}/ngmlr/sniffles/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/ngmlr/sniffles/total/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/sniffles_ngmlr_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/ngmlr/sniffles/total/tmp.txt"
    shell:
        """
        echo alignments/SI00001.ngmlr.srt.bam SNIFFLES > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10' |  bcftools view -e 'GT[*]="RR"' | bcftools sort > {output} && rm {params.tmp} 2> {log}
        """

rule clean_lra_sniffles:
    input:
        f"{RESULTDIR}/lra/sniffles/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/lra/sniffles/total/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/sniffles_lra_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/lra/sniffles/total/tmp.txt"
    shell:
        """
        echo alignments/SI00001.lra.srt.bam SNIFFLES > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10' |  bcftools view -e 'GT[*]="RR"' | bcftools sort > {output} && rm {params.tmp} 2> {log}
        """

rule clean_minimap2_svim:
    input:
        f"{RESULTDIR}/minimap2/svim/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/minimap2/svim/total/SI00001.clean.vcf"        
    log:
        f"{LOGDIR}/results/svim_minimap2_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/minimap2/svim/total/tmp.txt"
    shell:
        """
        echo SI00001 SVIM > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'SUPPORT>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output} && rm {params.tmp} 2> {log}
        """

rule clean_ngmlr_svim:
    input:
        f"{RESULTDIR}/ngmlr/svim/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/ngmlr/svim/total/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/svim_ngmlr_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/ngmlr/svim/total/tmp.txt"
    shell:
        """
        echo SI00001 SVIM > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'SUPPORT>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output} && rm {params.tmp} 2> {log}
        """

rule clean_lra_svim:
    input:
        f"{RESULTDIR}/lra/svim/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/lra/svim/total/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/svim_lra_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/lra/svim/total/tmp.txt"
    shell:
        """
        echo SI00001 SVIM > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'SUPPORT>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output} && rm {params.tmp} 2> {log}
        """

rule clean_pbmm2_pbsv:
    input:
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf"  
    log:
        f"{LOGDIR}/results/pbsv_pbmm2_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/pbmm2/pbsv/total/tmp.txt"
    shell:
        """
        echo SI00001 PBSV > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'AD[:1]>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output} && rm {params.tmp} 2> {log}
        """

rule clean_truth:
    input:
        f"{VCFTRUTH}"
    output:
        f"{VCFDIR}/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/truth_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{VCFDIR}/tmp.txt"
    shell:
        """
        echo SI00001 TRUTH > {params.tmp} && bcftools reheader -s {params.tmp} {input} | bcftools view -f PASS | grep -v -f {params.exclude}  | bcftools view -e 'GT[*]="RR"' | bcftools sort > {output} && rm {params.tmp} 2> {log}
        """

rule clean_minimap2_npinv:
    input:
        f"{RESULTDIR}/minimap2/npinv/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/minimap2/npinv/total/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/npinv_minimap2_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmpfile=f"{RESULTDIR}/minimap2/npinv/total/tmp.txt",
        tmpvcf=f"{RESULTDIR}/minimap2/npinv/total/tmp.vcf",
        tmpbed=f"{RESULTDIR}/minimap2/npinv/total/tmp.bed"
    shell:
        """
        echo Sample NPINV > {params.tmpfile} && grep -v -f {params.exclude} {input} | bcftools view -f PASS -i 'DV>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.tmpvcf} && bcftools query -f '%CHROM\t%POS\t%INFO/END\t%CHROM\t%POS\t%INFO/END\t%ID\t,\t+\t+\tINV\n' {params.tmpvcf} > {params.tmpbed} && rm {params.tmpvcf} && SURVIVOR bedtovcf {params.tmpbed} INV {params.tmpvcf} && rm {params.tmpbed} && bcftools reheader -s {params.tmpfile} {params.tmpvcf} > {output} && rm {params.tmpfile} && rm {params.tmpvcf} 2> {log}
        """

rule clean_ngmlr_npinv:
    input:
        f"{RESULTDIR}/ngmlr/npinv/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/ngmlr/npinv/total/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/npinv_ngmlr_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmpfile=f"{RESULTDIR}/ngmlr/npinv/total/tmp.txt",
        tmpvcf=f"{RESULTDIR}/ngmlr/npinv/total/tmp.vcf",
        tmpbed=f"{RESULTDIR}/ngmlr/npinv/total/tmp.bed"
    shell:
        """
        echo Sample NPINV > {params.tmpfile} && grep -v -f {params.exclude} {input} | bcftools view -f PASS -i 'DV>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.tmpvcf} && bcftools query -f '%CHROM\t%POS\t%INFO/END\t%CHROM\t%POS\t%INFO/END\t%ID\t,\t+\t+\tINV\n' {params.tmpvcf} > {params.tmpbed} && rm {params.tmpvcf} && SURVIVOR bedtovcf {params.tmpbed} INV {params.tmpvcf} && rm {params.tmpbed} && bcftools reheader -s {params.tmpfile} {params.tmpvcf} > {output} && rm {params.tmpfile} && rm {params.tmpvcf} 2> {log}
        """

rule clean_lra_npinv:
    input:
        f"{RESULTDIR}/lra/npinv/total/SI00001.vcf"
    output:
        f"{RESULTDIR}/lra/npinv/total/SI00001.clean.vcf"
    log:
        f"{LOGDIR}/results/npinv_lra_filter.log"
    conda: "../envs/sumsvs.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmpfile=f"{RESULTDIR}/lra/npinv/total/tmp.txt",
        tmpvcf=f"{RESULTDIR}/lra/npinv/total/tmp.vcf",
        tmpbed=f"{RESULTDIR}/lra/npinv/total/tmp.bed"
    shell:
        """
        echo Sample NPINV > {params.tmpfile} && grep -v -f {params.exclude} {input} | bcftools view -f PASS -i 'DV>=10' | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.tmpvcf} && bcftools query -f '%CHROM\t%POS\t%INFO/END\t%CHROM\t%POS\t%INFO/END\t%ID\t,\t+\t+\tINV\n' {params.tmpvcf} > {params.tmpbed} && rm {params.tmpvcf} && SURVIVOR bedtovcf {params.tmpbed} INV {params.tmpvcf} && rm {params.tmpbed} && bcftools reheader -s {params.tmpfile} {params.tmpvcf} > {output} && rm {params.tmpfile} && rm {params.tmpvcf} 2> {log}
        """


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
        SURVIVOR merge {params.tmp} 500 1 1 1 0 50 {params.mergedvcf} && rm {params.tmp}  && \
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
        SURVIVOR merge {params.tmp} 500 1 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
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
        SURVIVOR merge {params.tmp} 500 1 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
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
        SURVIVOR merge {params.tmp} 500 3 1 1 0 50 {output} && rm {params.tmp} 2> {log}
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
        SURVIVOR merge {params.tmp} 500 2 1 1 0 50 {output} && rm {params.tmp} 2> {log}
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
        SURVIVOR merge {params.tmp} 500 2 1 1 0 50 {output} && rm {params.tmp} 2> {log}
        """

rule survivor_consensus_all:
    input:
        expand(f"{RESULTDIR}/SI00001.{{aligner}}.consensus.vcf", aligner=["minimap2", "ngmlr", "lra"])
    output:
        f"{RESULTDIR}/SI00001.consensus.vcf"
    log:
        f"{LOGDIR}/results/survivor_consensus_all.log"
    threads: 1
    conda: "../envs/sumsvs.yaml"
    params:
        tmp=f"{RESULTDIR}/tmp_.txt"
    shell:
        """
        ls {input} > {params.tmp} && \
        SURVIVOR merge {params.tmp} 500 2 1 1 0 50 {output} && rm {params.tmp} 2> {log}
        """

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
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "cuteSV", "minimap2"}}' >> {output} 2>{log}
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
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "cuteSV", "NGMLR"}}' >> {output} 2>{log}
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
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "cuteSV", "lra"}}' >> {output} 2>{log}
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
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "Sniffles", "minimap2"}}' >> {output} 2>{log}
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
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "Sniffles", "NGMLR"}}' >> {output} 2>{log}
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
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "Sniffles", "lra"}}' >> {output} 2>{log}
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
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "pbsv", "minimap2"}}' >> {output} 2>{log}
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
        sed 's/DUP:TANDEM/DUP/g' {input} | sed 's/DUP:INT/DUP/g' > {params.tmp} && \
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
        sed 's/DUP:TANDEM/DUP/g' {input} | sed 's/DUP:INT/DUP/g' > {params.tmp} && \
        truvari stats -o {params.outstats} {params.tmp} && echo -e "SIZE\tDEL\tINS\tDUP\tINV\tTRA\ttool\taligner" > {output} && rm {params.tmp} && \
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "SVIM", "NGMLR"}}' >> {output} 2>{log}
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
        sed 's/DUP:TANDEM/DUP/g' {input} | sed 's/DUP:INT/DUP/g' > {params.tmp} && \
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
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "npInv", "minimap2"}}' >> {output} 2>{log}
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
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "npInv", "NGMLR"}}' >> {output} 2>{log}
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
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "npInv", "lra"}}' >> {output} 2>{log}
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
        cat {params.outstats} | grep -A 12 "# SVxSZ counts" | head -13 | tail -n +3 | sed 's/\[//g' | sed 's/)//g' | sed 's/,/-/g' | awk '{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$5,$7, "truth", "truth"}}' >> {output} 2>{log}
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

