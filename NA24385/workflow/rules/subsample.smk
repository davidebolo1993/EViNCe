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
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

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
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

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
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

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
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_30X:
    input:
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/30X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/30X/GM24385.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.6666667" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_minimap2_40X:
    input:
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/40X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/40X/GM24385.minimap2.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_subsample.log"
    conda: "../envs/minimap2.yaml"
    threads: 10
    params:
        fraction="0.8888889" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

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
        fraction="0.1111111" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

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
        fraction="0.2222222" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

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
        fraction="0.3333333" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

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
        fraction="0.4444444" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_30X:
    input:
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/30X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/30X/GM24385.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.6666667" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"

rule subsample_ngmlr_40X:
    input:
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        bam=f"{ALIGNDIR}/40X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/40X/GM24385.ngmlr.srt.bam.bai"
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_subsample.log"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    params:
        fraction="0.8888889" #considering 45X coverage as from mosdepth depth
    shell:
        "samtools view -s {params.fraction} -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam} 2>> {log}"
