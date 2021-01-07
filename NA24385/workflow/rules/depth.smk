rule mosdepth_minimap2_depth:
    input:
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        f"{ALIGNDIR}/GM24385.minimap2.mosdepth.global.dist.txt",
        f"{ALIGNDIR}/GM24385.minimap2.mosdepth.region.dist.txt",
        f"{ALIGNDIR}/GM24385.minimap2.mosdepth.summary.txt",
        f"{ALIGNDIR}/GM24385.minimap2.regions.bed.gz"
    log:
        f"{LOGDIR}/alignments/minimap2_mosdepth_depth.log"
    conda: "../envs/mosdepth.yaml"
    threads: 10
    params:
        prefix = f"{ALIGNDIR}/GM24385.minimap2"
    shell:
        "mosdepth -n --fast-mode -t {threads} {params.prefix} {input.bam} 2> {log}"

rule mosdepth_ngmlr_depth:
    input:
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{ALIGNDIR}/GM24385.ngmlr.mosdepth.global.dist.txt",
        f"{ALIGNDIR}/GM24385.ngmlr.mosdepth.region.dist.txt",
        f"{ALIGNDIR}/GM24385.ngmlr.mosdepth.summary.txt",
        f"{ALIGNDIR}/GM24385.ngmlr.regions.bed.gz"
    log:
        f"{LOGDIR}/alignments/ngmlr_mosdepth_depth.log"
    conda: "../envs/mosdepth.yaml"
    threads: 10
    params:
        prefix = f"{ALIGNDIR}/GM24385.ngmlr"
    shell:
        "mosdepth -n --fast-mode -t {threads} {params.prefix} {input.bam} 2> {log}"

rule plot_cov:
    input:
        expand(f"{ALIGNDIR}/GM24385.{{aligner}}.mosdepth.global.dist.txt", aligner=['minimap2','ngmlr']),
	expand(f"{ALIGNDIR}/GM24385.{{aligner}}.regions.bed.gz", aligner=['minimap2','ngmlr'])
    output:
        expand(f"{ALIGNDIR}/GM24385.{{aligner}}.coverage_per_chromosome.png", aligner=['minimap2','ngmlr']),
	f"{ALIGNDIR}/GM24385.fraction_of_bases_per_depth.png"
    log:
        f"{LOGDIR}/alignments/plot_depth.log"
    conda:
        "../envs/plot.yaml"
    params:
        indir=f"{ALIGNDIR}",
        script=f"../scripts/coveragestats.R"
    shell:
        "Rscript {params.script} {params.indir} 2> {log}" 