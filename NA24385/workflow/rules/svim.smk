rule svim_call_total_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/total/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_total_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/total",
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}" 

rule svim_call_total_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/total/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_total_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/total"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}"

rule svim_call_5X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/5X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_5X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/5X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}" 

rule svim_call_5X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/5X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_5X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/5X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}"

rule svim_call_10X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/10X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_10X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/10X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}" 

rule svim_call_10X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/10X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_10X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/10X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}"

rule svim_call_15X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/15X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_15X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/15X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}" 

rule svim_call_15X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/15X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_15X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/15X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}"

rule svim_call_20X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/20X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_20X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/20X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}" 

rule svim_call_20X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/20X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_20X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/20X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}"


rule svim_call_25X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/25X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_25X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/25X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}" 

rule svim_call_25X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/25X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_25X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/25X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}"


rule svim_call_35X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/svim/35X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_svim_35X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/svim/35X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}" 

rule svim_call_35X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/svim/35X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_svim_35X_call.log"
    conda: "../envs/svim.yaml"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/svim/35X"
    shell:
        "mkdir -p {params.wd} && svim alignment --min_sv_size 50 --sample GM24385 {params.wd} {input.bam} {input.ref} && mv {params.wd}/varants.vcf {params.wd}/GM24385.vcf 2>{log}"
