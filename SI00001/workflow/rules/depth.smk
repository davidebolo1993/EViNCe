rule mosdepth_minimap2_depth:
    input:
        bam=f"{ALIGNDIR}/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    output:
        f"{ALIGNDIR}/SI00001.minimap2.mosdepth.global.dist.txt",
        f"{ALIGNDIR}/SI00001.minimap2.mosdepth.region.dist.txt",
        f"{ALIGNDIR}/SI00001.minimap2.mosdepth.summary.txt",
        f"{ALIGNDIR}/SI00001.minimap2.regions.bed.gz",
        f"{ALIGNDIR}/SI00001.minimap2.regions.bed.gz.csi"
    log:
        f"{LOGDIR}/alignments/minimap2_mosdepth_depth.log"
    conda: "../envs/mosdepth.yaml"
    threads: 10
    params:
        prefix = f"{ALIGNDIR}/SI00001.minimap2"
    shell:
        "mosdepth -n --fast-mode -t {threads} --by 500 {params.prefix} {input.bam} 2> {log}"

rule mosdepth_ngmlr_depth:
    input:
        bam=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{ALIGNDIR}/SI00001.ngmlr.mosdepth.global.dist.txt",
        f"{ALIGNDIR}/SI00001.ngmlr.mosdepth.region.dist.txt",
        f"{ALIGNDIR}/SI00001.ngmlr.mosdepth.summary.txt",
        f"{ALIGNDIR}/SI00001.ngmlr.regions.bed.gz",
        f"{ALIGNDIR}/SI00001.ngmlr.regions.bed.gz.csi"
    log:
        f"{LOGDIR}/alignments/ngmlr_mosdepth_depth.log"
    conda: "../envs/mosdepth.yaml"
    threads: 10
    params:
        prefix = f"{ALIGNDIR}/SI00001.ngmlr"
    shell:
        "mosdepth -n --fast-mode -t {threads} --by 500 {params.prefix} {input.bam} 2> {log}"

rule mosdepth_lra_depth:
    input:
        bam=f"{ALIGNDIR}/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    output:
        f"{ALIGNDIR}/SI00001.lra.mosdepth.global.dist.txt",
        f"{ALIGNDIR}/SI00001.lra.mosdepth.region.dist.txt",
        f"{ALIGNDIR}/SI00001.lra.mosdepth.summary.txt",
        f"{ALIGNDIR}/SI00001.lra.regions.bed.gz",
        f"{ALIGNDIR}/SI00001.lra.regions.bed.gz.csi"
    log:
        f"{LOGDIR}/alignments/lra_mosdepth_depth.log"
    conda: "../envs/mosdepth.yaml"
    threads: 10
    params:
        prefix = f"{ALIGNDIR}/SI00001.lra"
    shell:
        "mosdepth -n --fast-mode -t {threads} --by 500 {params.prefix} {input.bam} 2> {log}"

rule mosdepth_pbmm2_depth:
    input:
        bam=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam.bai"
    output:
        f"{ALIGNDIR}/SI00001.pbmm2.mosdepth.global.dist.txt",
        f"{ALIGNDIR}/SI00001.pbmm2.mosdepth.region.dist.txt",
        f"{ALIGNDIR}/SI00001.pbmm2.mosdepth.summary.txt",
        f"{ALIGNDIR}/SI00001.pbmm2.regions.bed.gz",
        f"{ALIGNDIR}/SI00001.pbmm2.regions.bed.gz.csi"
    log:
        f"{LOGDIR}/alignments/pbmm2_mosdepth_depth.log"
    conda: "../envs/mosdepth.yaml"
    threads: 10
    params:
        prefix = f"{ALIGNDIR}/SI00001.pbmm2"
    shell:
        "mosdepth -n --fast-mode -t {threads} --by 500 {params.prefix} {input.bam} 2> {log}"

rule plot_cov:
    input:
        expand(f"{ALIGNDIR}/SI00001.{{aligner}}.mosdepth.global.dist.txt", aligner=['minimap2','pbmm2','ngmlr','lra'])
    output:
       f"{ALIGNDIR}/SI00001.fraction_of_bases_per_depth.pdf"
    log:
        f"{LOGDIR}/alignments/plot_depth.log"
    conda:
        "../envs/plot.yaml"
    shell:
        "Rscript {SCRIPTDIR}/coveragestats.R {ALIGNDIR} {VCFDIR} 2> {log}"
