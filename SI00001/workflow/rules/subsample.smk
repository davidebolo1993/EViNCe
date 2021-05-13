rule subsample_minimap2_5X:
    input:
        bam=f"{ALIGNDIR}/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/5X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.1" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_10X:
    input:
        bam=f"{ALIGNDIR}/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/10X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.2" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_15X:
    input:
        bam=f"{ALIGNDIR}/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/15X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.3" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_20X:
    input:
        bam=f"{ALIGNDIR}/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/20X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.4" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_25X:
    input:
        bam=f"{ALIGNDIR}/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/25X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.5" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_35X:
    input:
        bam=f"{ALIGNDIR}/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/35X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.7" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_5X:
    input:
        bam=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/5X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.1" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_10X:
    input:
        bam=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/10X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.2" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_15X:
    input:
        bam=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/15X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.3" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_20X:
    input:
        bam=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/20X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.4" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_25X:
    input:
        bam=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/25X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.5" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_35X:
    input:
        bam=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/35X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.7" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_5X:
    input:
        bam=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/5X/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.1" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_10X:
    input:
        bam=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/10X/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.2" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_15X:
    input:
        bam=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/15X/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.3" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_20X:
    input:
        bam=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/20X/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.4" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_25X:
    input:
        bam=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/25X/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.5" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_35X:
    input:
        bam=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/35X/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.7" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_5X:
    input:
        bam=f"{ALIGNDIR}/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/5X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.1" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_10X:
    input:
        bam=f"{ALIGNDIR}/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/10X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.2" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_15X:
    input:
        bam=f"{ALIGNDIR}/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/15X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.3" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_20X:
    input:
        bam=f"{ALIGNDIR}/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/20X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.4" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_25X:
    input:
        bam=f"{ALIGNDIR}/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/25X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.5" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_35X:
    input:
        bam=f"{ALIGNDIR}/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/35X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.7" #considering 50X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"
