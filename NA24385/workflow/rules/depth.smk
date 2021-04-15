rule mosdepth_minimap2_depth:
    input:
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        f"{ALIGNDIR}/GM24385.minimap2.mosdepth.global.dist.txt",
        f"{ALIGNDIR}/GM24385.minimap2.mosdepth.region.dist.txt",
        f"{ALIGNDIR}/GM24385.minimap2.mosdepth.summary.txt",
        f"{ALIGNDIR}/GM24385.minimap2.regions.bed.gz",
        f"{ALIGNDIR}/GM24385.minimap2.regions.bed.gz.csi"
    log:
        f"{LOGDIR}/alignments/minimap2_mosdepth_depth.log"
    conda: "../envs/mosdepth.yaml"
    threads: 10
    params:
        prefix = f"{ALIGNDIR}/GM24385.minimap2"
    shell:
        "mosdepth -n --fast-mode -t {threads} --by 500 {params.prefix} {input.bam} 2> {log}"

rule mosdepth_ngmlr_depth:
    input:
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{ALIGNDIR}/GM24385.ngmlr.mosdepth.global.dist.txt",
        f"{ALIGNDIR}/GM24385.ngmlr.mosdepth.region.dist.txt",
        f"{ALIGNDIR}/GM24385.ngmlr.mosdepth.summary.txt",
        f"{ALIGNDIR}/GM24385.ngmlr.regions.bed.gz",
        f"{ALIGNDIR}/GM24385.ngmlr.regions.bed.gz.csi"
    log:
        f"{LOGDIR}/alignments/ngmlr_mosdepth_depth.log"
    conda: "../envs/mosdepth.yaml"
    threads: 10
    params:
        prefix = f"{ALIGNDIR}/GM24385.ngmlr"
    shell:
        "mosdepth -n --fast-mode -t {threads} --by 500 {params.prefix} {input.bam} 2> {log}"

rule mosdepth_lra_depth:
    input:
        bam=f"{ALIGNDIR}/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.lra.srt.bam.bai"
    output:
        f"{ALIGNDIR}/GM24385.lra.mosdepth.global.dist.txt",
        f"{ALIGNDIR}/GM24385.lra.mosdepth.region.dist.txt",
        f"{ALIGNDIR}/GM24385.lra.mosdepth.summary.txt",
        f"{ALIGNDIR}/GM24385.lra.regions.bed.gz",
        f"{ALIGNDIR}/GM24385.lra.regions.bed.gz.csi"
    log:
        f"{LOGDIR}/alignments/lra_mosdepth_depth.log"
    conda: "../envs/mosdepth.yaml"
    threads: 10
    params:
        prefix = f"{ALIGNDIR}/GM24385.lra"
    shell:
        "mosdepth -n --fast-mode -t {threads} --by 500 {params.prefix} {input.bam} 2> {log}"

rule mosdepth_pbmm2_depth:
    input:
        bam=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam.bai"
    output:
        f"{ALIGNDIR}/GM24385.pbmm2.mosdepth.global.dist.txt",
        f"{ALIGNDIR}/GM24385.pbmm2.mosdepth.region.dist.txt",
        f"{ALIGNDIR}/GM24385.pbmm2.mosdepth.summary.txt",
        f"{ALIGNDIR}/GM24385.pbmm2.regions.bed.gz",
        f"{ALIGNDIR}/GM24385.pbmm2.regions.bed.gz.csi"
    log:
        f"{LOGDIR}/alignments/pbmm2_mosdepth_depth.log"
    conda: "../envs/mosdepth.yaml"
    threads: 10
    params:
        prefix = f"{ALIGNDIR}/GM24385.pbmm2"
    shell:
        "mosdepth -n --fast-mode -t {threads} --by 500 {params.prefix} {input.bam} 2> {log}"

rule plot_cov:
    input:
        expand(f"{ALIGNDIR}/GM24385.{{aligner}}.mosdepth.global.dist.txt", aligner=['minimap2','pbmm2','ngmlr','lra']),
        expand(f"{ALIGNDIR}/GM24385.{{aligner}}.regions.bed.gz", aligner=['minimap2','pbmm2','ngmlr','lra'])
    output:
        expand(f"{ALIGNDIR}/GM24385.{{aligner}}.coverage_per_chromosome.png", aligner=['minimap2','pbmm2','ngmlr','lra']),
        expand(f"{VCFDIR}/GM24385.{{aligner}}.exclude.tsv", aligner=['minimap2','pbmm2','ngmlr','lra']),
        f"{ALIGNDIR}/GM24385.fraction_of_bases_per_depth.png"
    log:
        f"{LOGDIR}/alignments/plot_depth.log"
    conda:
        "../envs/plot.yaml"
    shell:
        "Rscript {SCRIPTDIR}/coveragestats.R {ALIGNDIR} {VCFDIR} 2> {log}"


