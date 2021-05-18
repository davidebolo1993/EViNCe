rule sniffles_call_total_minimap2:
    input:
        bam=f"{ALIGNDIR}/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/sniffles/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_sniffles_total_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/minimap2/sniffles/total/sniffles.tmp",
        wd=f"{RESULTDIR}/minimap2/sniffles/total"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}" 

rule sniffles_call_total_ngmlr:
    input:
        bam=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/sniffles/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_sniffles_total_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/ngmlr/sniffles/total/sniffles.tmp",
        wd=f"{RESULTDIR}/ngmlr/sniffles/total"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_total_lra:
    input:
        bam=f"{ALIGNDIR}/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/sniffles/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_sniffles_total_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/lra/sniffles/total/sniffles.tmp",
        wd=f"{RESULTDIR}/lra/sniffles/total"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_5X_minimap2:
    input:
        bam=f"{ALIGNDIR}/5X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/sniffles/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_sniffles_5X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/minimap2/sniffles/5X/sniffles.tmp",
        wd=f"{RESULTDIR}/minimap2/sniffles/5X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}" 

rule sniffles_call_5X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/5X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/sniffles/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_sniffles_5X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/ngmlr/sniffles/5X/sniffles.tmp",
        wd=f"{RESULTDIR}/ngmlr/sniffles/5X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_5X_lra:
    input:
        bam=f"{ALIGNDIR}/5X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/sniffles/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_sniffles_5X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/lra/sniffles/5X/sniffles.tmp",
        wd=f"{RESULTDIR}/lra/sniffles/5X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_10X_minimap2:
    input:
        bam=f"{ALIGNDIR}/10X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/sniffles/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_sniffles_10X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/minimap2/sniffles/10X/sniffles.tmp",
        wd=f"{RESULTDIR}/minimap2/sniffles/10X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}" 

rule sniffles_call_10X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/10X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/sniffles/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_sniffles_10X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/ngmlr/sniffles/10X/sniffles.tmp",
        wd=f"{RESULTDIR}/ngmlr/sniffles/10X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_10X_lra:
    input:
        bam=f"{ALIGNDIR}/10X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/sniffles/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_sniffles_10X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/lra/sniffles/10X/sniffles.tmp",
        wd=f"{RESULTDIR}/lra/sniffles/10X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"


rule sniffles_call_15X_minimap2:
    input:
        bam=f"{ALIGNDIR}/15X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/sniffles/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_sniffles_15X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/minimap2/sniffles/15X/sniffles.tmp",
        wd=f"{RESULTDIR}/minimap2/sniffles/15X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}" 

rule sniffles_call_15X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/15X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/sniffles/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_sniffles_15X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/ngmlr/sniffles/15X/sniffles.tmp",
        wd=f"{RESULTDIR}/ngmlr/sniffles/15X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_15X_lra:
    input:
        bam=f"{ALIGNDIR}/15X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/sniffles/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_sniffles_15X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/lra/sniffles/15X/sniffles.tmp",
        wd=f"{RESULTDIR}/lra/sniffles/15X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_20X_minimap2:
    input:
        bam=f"{ALIGNDIR}/20X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/sniffles/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_sniffles_20X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/minimap2/sniffles/20X/sniffles.tmp",
        wd=f"{RESULTDIR}/minimap2/sniffles/20X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}" 

rule sniffles_call_20X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/20X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/sniffles/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_sniffles_20X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/ngmlr/sniffles/20X/sniffles.tmp",
        wd=f"{RESULTDIR}/ngmlr/sniffles/20X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_20X_lra:
    input:
        bam=f"{ALIGNDIR}/20X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/sniffles/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_sniffles_20X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/lra/sniffles/20X/sniffles.tmp",
        wd=f"{RESULTDIR}/lra/sniffles/20X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_25X_minimap2:
    input:
        bam=f"{ALIGNDIR}/25X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/sniffles/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_sniffles_25X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/minimap2/sniffles/25X/sniffles.tmp",
        wd=f"{RESULTDIR}/minimap2/sniffles/25X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}" 

rule sniffles_call_25X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/25X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/sniffles/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_sniffles_25X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/ngmlr/sniffles/25X/sniffles.tmp",
        wd=f"{RESULTDIR}/ngmlr/sniffles/25X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_25X_lra:
    input:
        bam=f"{ALIGNDIR}/25X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/sniffles/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_sniffles_25X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/lra/sniffles/25X/sniffles.tmp",
        wd=f"{RESULTDIR}/lra/sniffles/25X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_35X_minimap2:
    input:
        bam=f"{ALIGNDIR}/35X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/sniffles/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_sniffles_35X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/minimap2/sniffles/35X/sniffles.tmp",
        wd=f"{RESULTDIR}/minimap2/sniffles/35X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}" 

rule sniffles_call_35X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/35X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/sniffles/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_sniffles_35X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/ngmlr/sniffles/35X/sniffles.tmp",
        wd=f"{RESULTDIR}/ngmlr/sniffles/35X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

rule sniffles_call_35X_lra:
    input:
        bam=f"{ALIGNDIR}/35X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/sniffles/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_sniffles_35X_call.log"
    conda: "../envs/sniffles.yaml"
    threads: 10
    params:
        wf=f"{RESULTDIR}/lra/sniffles/35X/sniffles.tmp",
        wd=f"{RESULTDIR}/lra/sniffles/35X"
    shell:
        "mkdir -p {params.wd} && sniffles -m {input.bam} -v {output} --tmp_file {params.wf} -s 2 -l 50 -t {threads} 2>{log}"

