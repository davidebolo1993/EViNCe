rule clean_minimap2_cutesv:
    input:
        f"{RESULTDIR}/minimap2/cutesv/total/GM24385.vcf",
    output:
        f"{RESULTDIR}/minimap2/cutesv/total/GM24385.clean.vcf",
    log:
        f"{LOGDIR}/results/cutesv_minimap2_filter.log"
    conda: "../envs/upset.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/minimap2/cutesv/total/tmp.txt"
    shell:
        "echo GM24385 CUTESV > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10'| bcftools sort > {output} && rm {params.tmp} 2> {log}"

rule clean_ngmlr_cutesv:
    input:
        f"{RESULTDIR}/ngmlr/cutesv/total/GM24385.vcf",
    output:
        f"{RESULTDIR}/ngmlr/cutesv/total/GM24385.clean.vcf",
    log:
        f"{LOGDIR}/results/cutesv_ngmlr_filter.log"
    conda: "../envs/upset.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/ngmlr/cutesv/total/tmp.txt"
    shell:
        "echo GM24385 CUTESV > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10'| bcftools sort > {output} && rm {params.tmp} 2> {log}"

rule clean_minimap2_sniffles:
    input:
        f"{RESULTDIR}/minimap2/sniffles/total/GM24385.vcf",
    output:
        f"{RESULTDIR}/minimap2/sniffles/total/GM24385.clean.vcf",
    log:
        f"{LOGDIR}/results/sniffles_minimap2_filter.log"
    conda: "../envs/upset.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/minimap2/sniffles/total/tmp.txt"
    shell:
        "echo alignments/GM24385.minimap2.srt.bam SNIFFLES > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10' | bcftools sort > {output} && rm {params.tmp} 2> {log}"

rule clean_ngmlr_sniffles:
    input:
        f"{RESULTDIR}/ngmlr/sniffles/total/GM24385.vcf",
    output:
        f"{RESULTDIR}/ngmlr/sniffles/total/GM24385.clean.vcf",
    log:
        f"{LOGDIR}/results/sniffles_ngmlr_filter.log"
    conda: "../envs/upset.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/ngmlr/sniffles/total/tmp.txt"
    shell:
        "echo alignments/GM24385.ngmlr.srt.bam SNIFFLES > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'DV>=10' | bcftools sort > {output} && rm {params.tmp} 2> {log}"

rule clean_minimap2_svim:
    input:
        f"{RESULTDIR}/minimap2/svim/total/GM24385.vcf",
    output:
        f"{RESULTDIR}/minimap2/svim/total/GM24385.clean.vcf",
    log:
        f"{LOGDIR}/results/svim_minimap2_filter.log"
    conda: "../envs/upset.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/minimap2/svim/total/tmp.txt"
    shell:
        "echo GM24385 SVIM > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'SUPPORT>=10'| bcftools sort > {output} && rm {params.tmp} 2> {log}"

rule clean_ngmlr_svim:
    input:
        f"{RESULTDIR}/ngmlr/svim/total/GM24385.vcf",
    output:
        f"{RESULTDIR}/ngmlr/svim/total/GM24385.clean.vcf",
    log:
        f"{LOGDIR}/results/svim_ngmlr_filter.log"
    conda: "../envs/upset.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/ngmlr/svim/total/tmp.txt"
    shell:
        "echo GM24385 SVIM > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'SUPPORT>=10' | bcftools sort > {output} && rm {params.tmp} 2> {log}"

rule clean_minimap2_npinv:
    input:
        f"{RESULTDIR}/minimap2/npinv/total/GM24385.vcf",
    output:
        f"{RESULTDIR}/minimap2/npinv/total/GM24385.clean.vcf",
    log:
        f"{LOGDIR}/results/npinv_minimap2_filter.log"
    conda: "../envs/upset.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmpfile=f"{RESULTDIR}/minimap2/npinv/total/tmp.txt",
        tmpvcf=f"{RESULTDIR}/minimap2/npinv/total/tmp.vcf",
        tmpbed=f"{RESULTDIR}/minimap2/npinv/total/tmp.bed",
    shell:
        "echo Sample NPINV > {params.tmpfile} && grep -v -f {params.exclude} {input} | bcftools view -f PASS -i 'DV>=10'| bcftools sort > {params.tmpvcf} && bcftools query -f '%CHROM\t%POS\t%INFO/END\t%CHROM\t%POS\t%INFO/END\t%ID\t,\t+\t+\tINV\n' {params.tmpvcf} > {params.tmpbed} && rm {params.tmpvcf} && SURVIVOR bedtovcf {params.tmpbed} INV {params.tmpvcf} && rm {params.tmpbed} && bcftools reheader -s {params.tmpfile} {params.tmpvcf} > {output} && rm {params.tmpfile} 2> {log}"

rule clean_ngmlr_npinv:
    input:
        f"{RESULTDIR}/ngmlr/npinv/total/GM24385.vcf",
    output:
        f"{RESULTDIR}/ngmlr/npinv/total/GM24385.clean.vcf",
    log:
        f"{LOGDIR}/results/npinv_ngmlr_filter.log"
    conda: "../envs/upset.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmpfile=f"{RESULTDIR}/ngmlr/npinv/total/tmp.txt",
        tmpvcf=f"{RESULTDIR}/ngmlr/npinv/total/tmp.vcf",
        tmpbed=f"{RESULTDIR}/ngmlr/npinv/total/tmp.bed",
    shell:
        "echo Sample NPINV > {params.tmpfile} && grep -v -f {params.exclude} {input} | bcftools view -f PASS -i 'DV>=10' | bcftools sort > {params.tmpvcf} && bcftools query -f '%CHROM\t%POS\t%INFO/END\t%CHROM\t%POS\t%INFO/END\t%ID\t,\t+\t+\tINV\n' {params.tmpvcf} > {params.tmpbed} && rm {params.tmpvcf} && SURVIVOR bedtovcf {params.tmpbed} INV {params.tmpvcf} && rm {params.tmpbed} && bcftools reheader -s {params.tmpfile} {params.tmpvcf} > {output} && rm {params.tmpfile} 2> {log}"

rule clean_pbmm2_pbsv:
    input:
        f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.vcf",
    output:
        f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.clean.vcf",
    log:
        f"{LOGDIR}/results/pbsv_pbmm2_filter.log"
    conda: "../envs/upset.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{RESULTDIR}/pbmm2/pbsv/total/tmp.txt"
    shell:
        "echo GM24385 PBSV > {params.tmp} && bcftools reheader -s {params.tmp} {input} | grep -v -f {params.exclude} | bcftools view -f PASS -i 'SAC[:1]>=10' | bcftools sort > {output} && rm {params.tmp} 2> {log}"

rule clean_truth:
    input:
        f"{VCFTRUTH}",
    output:
        f"{VCFDIR}/GM24385.clean.vcf",
    log:
        f"{LOGDIR}/results/truth_filter.log"
    conda: "../envs/upset.yaml"
    threads: 1
    params:
        exclude = config["excludebed"],
        tmp=f"{VCFDIR}/tmp.txt"
    shell:
        "echo HG002 TRUTH > {params.tmp} && bcftools reheader -s {params.tmp} {input} | bcftools view -f PASS | grep -v -f {params.exclude}  | bcftools sort > {output} && rm {params.tmp} 2> {log}"

rule survivor_merge_minimap2:
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv", "sniffles", "svim", "npinv"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.clean.vcf",
        f"{VCFDIR}/GM24385.clean.vcf"
    output:
        f"{RESULTDIR}/GM24385.minimap2.upset.tsv"
    log:
        f"{LOGDIR}/results/minimap2_survivor.log"
    threads: 1
    conda: "../envs/upset.yaml"
    params:
        tmp=f"{RESULTDIR}/tmpminimap2.txt",
        mergedvcf=f"{RESULTDIR}/GM24385.minimap2.merged.vcf"
    shell:
        """
        ls {input} > {params.tmp} && \
        SURVIVOR merge {params.tmp} 1000 1 1 1 0 50 {params.mergedvcf} && rm {params.tmp}  && \
        bcftools view -h {params.mergedvcf} | tail -1 | awk '{{OFS=FS="\t"}}{{print $3,$10,$11,$12,$13,$14,$15,$5,$1}}' > {output}  && \
        bcftools query -f '%CHROM\t%SVTYPE[\t%GT:%PSV:%LN:%DR:%ST:%QV:%TY:%ID:%RAL:%AAL:%CO]\n' {params.mergedvcf}  | awk '{{FS=OFS="\t"}}{{print "var"NR,$3,$4,$5,$6,$7,$8,$2,$1}}' | awk '{{FS=OFS="\t"}}{{if ($2 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $2="0"; else $2="1"; if ($3 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $3="0"; else $3="1"; if ($4 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $4="0"; else $4="1"; if ($5 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $5="0"; else $5="1";if ($6 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $6="0"; else $6="1"; if ($7 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $7="0"; else $7="1"}}1' >> {output} 2>{log}
        """

rule survivor_merge_ngmlr:
    input:
        expand(f"{RESULTDIR}/ngmlr/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv", "sniffles", "svim", "npinv"]),
        f"{VCFDIR}/GM24385.clean.vcf"
    output:
        f"{RESULTDIR}/GM24385.ngmlr.upset.tsv"
    log:
        f"{LOGDIR}/results/ngmlr_survivor.log"
    threads: 1
    conda: "../envs/upset.yaml"
    params:
        tmp=f"{RESULTDIR}/tmpngmlr.txt",
        mergedvcf=f"{RESULTDIR}/GM24385.ngmlr.merged.vcf"
    shell:
        """
        ls {input} > {params.tmp} && \
        SURVIVOR merge {params.tmp} 1000 1 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -h {params.mergedvcf} | tail -1 |  awk '{{OFS=FS="\t"}}{{print $3,$10,$11,$12,$13,$14,$5,$1}}' > {output}  && \
        bcftools query -f '%CHROM\t%SVTYPE[\t%GT:%PSV:%LN:%DR:%ST:%QV:%TY:%ID:%RAL:%AAL:%CO]\n' {params.mergedvcf} | awk '{{FS=OFS="\t"}}{{print "var"NR,$3,$4,$5,$6,$7,$2,$1}}' | awk '{{FS=OFS="\t"}}{{if ($2 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $2="0"; else $2="1"; if ($3 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $3="0"; else $3="1"; if ($4 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $4="0"; else $4="1"; if ($5 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $5="0"; else $5="1";if ($6 == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN") $6="0"; else $6="1"}}1' >> {output} 2>{log}
        """

rule upset_plot:
    input:
        expand(f"{RESULTDIR}/GM24385.{{aligner}}.upset.tsv", aligner=["minimap2", "ngmlr"])
    output:
        expand(f"{RESULTDIR}/GM24385.{{aligner}}.upset.pdf", aligner=["minimap2", "ngmlr"])
    log:
        f"{LOGDIR}/results/upsetplot.log"
    threads: 1
    conda: "../envs/plot.yaml"
    shell:
        "Rscript {SCRIPTDIR}/upsetplot.R {RESULTDIR} 2> {log}"
