rule svim_call_total_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_total_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/total",
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}" 

rule svim_call_total_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_total_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/total"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"

rule svim_call_total_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/svim/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_svim_total_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/svim/total"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"

rule svim_call_5X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_5X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/5X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}" 

rule svim_call_5X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_5X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/5X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"

rule svim_call_5X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/svim/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_svim_5X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/svim/5X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"

rule svim_call_10X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_10X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/10X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}" 

rule svim_call_10X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_10X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/10X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"

rule svim_call_10X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/svim/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_svim_10X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/svim/10X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"


rule svim_call_15X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_15X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/15X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}" 

rule svim_call_15X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_15X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/15X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"


rule svim_call_15X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/svim/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_svim_15X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/svim/15X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"


rule svim_call_20X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_20X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/20X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}" 

rule svim_call_20X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_20X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/20X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"

rule svim_call_20X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/svim/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_svim_20X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/svim/20X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"


rule svim_call_25X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_25X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/25X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}" 

rule svim_call_25X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_25X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/25X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"

rule svim_call_25X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/svim/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_svim_25X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/svim/25X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"


rule svim_call_35X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_35X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/35X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}" 

rule svim_call_35X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_35X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/35X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"

rule svim_call_35X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/svim/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_svim_35X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/lra/svim/35X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample SI00001 {params.wd} {input.bam} {input.ref} && mv {params.wd}/variants.vcf {params.wd}/SI00001.vcf 2>{log}"


