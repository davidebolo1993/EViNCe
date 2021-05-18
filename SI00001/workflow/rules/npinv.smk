rule npinv_call_total_minimap2:
    input:
        bam=f"{ALIGNDIR}/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/npinv/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_npinv_total_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/npinv/total",
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}" 

rule npinv_call_total_ngmlr:
    input:
        bam=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/npinv/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_npinv_total_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/npinv/total"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"

rule npinv_call_total_lra:
    input:
        bam=f"{ALIGNDIR}/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/npinv/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_npinv_total_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/npinv/total"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"

rule npinv_call_5X_minimap2:
    input:
        bam=f"{ALIGNDIR}/5X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/npinv/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_npinv_5X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/npinv/5X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}" 

rule npinv_call_5X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/5X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/npinv/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_npinv_5X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/npinv/5X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"

rule npinv_call_5X_lra:
    input:
        bam=f"{ALIGNDIR}/5X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/npinv/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_npinv_5X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/npinv/5X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"

rule npinv_call_10X_minimap2:
    input:
        bam=f"{ALIGNDIR}/10X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/npinv/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_npinv_10X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/npinv/10X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}" 

rule npinv_call_10X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/10X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/npinv/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_npinv_10X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/npinv/10X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"

rule npinv_call_10X_lra:
    input:
        bam=f"{ALIGNDIR}/10X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/npinv/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_npinv_10X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/npinv/10X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"


rule npinv_call_15X_minimap2:
    input:
        bam=f"{ALIGNDIR}/15X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/npinv/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_npinv_15X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/npinv/15X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}" 

rule npinv_call_15X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/15X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/npinv/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_npinv_15X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/npinv/15X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"

rule npinv_call_15X_lra:
    input:
        bam=f"{ALIGNDIR}/15X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/npinv/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_npinv_15X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/npinv/15X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"


rule npinv_call_20X_minimap2:
    input:
        bam=f"{ALIGNDIR}/20X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/npinv/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_npinv_20X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/npinv/20X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}" 

rule npinv_call_20X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/20X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/npinv/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_npinv_20X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/npinv/20X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"

rule npinv_call_20X_lra:
    input:
        bam=f"{ALIGNDIR}/20X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/npinv/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_npinv_20X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/npinv/20X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"

rule npinv_call_25X_minimap2:
    input:
        bam=f"{ALIGNDIR}/25X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/npinv/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_npinv_25X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/npinv/25X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}" 

rule npinv_call_25X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/25X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/npinv/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_npinv_25X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/npinv/25X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"

rule npinv_call_25X_lra:
    input:
        bam=f"{ALIGNDIR}/25X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/npinv/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_npinv_25X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/npinv/25X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"

rule npinv_call_35X_minimap2:
    input:
        bam=f"{ALIGNDIR}/35X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/npinv/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_npinv_35X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/npinv/35X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}" 

rule npinv_call_35X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/35X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/npinv/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_npinv_35X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/npinv/35X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"

rule npinv_call_35X_lra:
    input:
        bam=f"{ALIGNDIR}/35X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/npinv/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_npinv_35X_call.log"
    conda: "../envs/npinv.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/npinv/35X"
    shell:
        "mkdir -p {params.wd} && npinv --output {output} --input {input.bam} --min 50 --threshold 2 2>{log}"


