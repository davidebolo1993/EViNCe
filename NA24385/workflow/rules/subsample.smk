rule subsample_minimap2_5X:
    input:
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/5X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.1111111" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_10X:
    input:
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/10X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.2222222" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_15X:
    input:
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/15X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.3333333" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_20X:
    input:
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/20X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.4444444" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_25X:
    input:
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/25X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.5555556" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_35X:
    input:
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/35X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.7777778" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_5X:
    input:
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/5X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.1190476" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_10X:
    input:
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/10X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.2380952" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_15X:
    input:
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/15X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.3571429" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_20X:
    input:
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/20X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.4761905" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_25X:
    input:
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/25X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.5952381" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_35X:
    input:
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/35X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.8333333" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_5X:
    input:
        bam=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/5X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.1190476" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_10X:
    input:
        bam=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/10X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.2380952" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_15X:
    input:
        bam=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/15X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.3571429" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_20X:
    input:
        bam=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/20X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.4761905" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_25X:
    input:
        bam=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/25X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.5952381" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_pbmm2_35X:
    input:
        bam=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/35X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.pbmm2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/pbmm2_samtools_subsample.log"
    conda: "../envs/pbmm2.yaml"
    threads: 10
    params:
        fraction="0.8333333" #considering 42X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_5X:
    input:
        bam=f"{ALIGNDIR}/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/5X/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.125" #considering 40X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_10X:
    input:
        bam=f"{ALIGNDIR}/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/10X/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.25" #considering 40X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_15X:
    input:
        bam=f"{ALIGNDIR}/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/15X/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.375" #considering 40X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_20X:
    input:
        bam=f"{ALIGNDIR}/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/20X/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.5" #considering 40X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_25X:
    input:
        bam=f"{ALIGNDIR}/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/25X/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.625" #considering 40X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_lra_35X:
    input:
        bam=f"{ALIGNDIR}/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.lra.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/35X/GM24385.lra.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.lra.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/lra_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.875" #considering 40X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"
