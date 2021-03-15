rule nanomonsv_call_total_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanomonsv/total/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanomonsv_total_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/nanomonsv/total",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_total_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanomonsv/total/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanomonsv_total_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/nanomonsv/total",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_35X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanomonsv/35X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanomonsv_35X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/nanomonsv/35X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_35X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanomonsv/35X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanomonsv_35X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/nanomonsv/35X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_25X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanomonsv/25X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanomonsv_25X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/nanomonsv/25X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_25X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanomonsv/25X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanomonsv_25X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/nanomonsv/25X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"


rule nanomonsv_call_20X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanomonsv/20X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanomonsv_20X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/nanomonsv/20X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_20X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanomonsv/20X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanomonsv_20X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/nanomonsv/20X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_15X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanomonsv/15X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanomonsv_15X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/nanomonsv/15X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_15X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanomonsv/15X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanomonsv_15X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/nanomonsv/15X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_10X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanomonsv/10X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanomonsv_10X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/nanomonsv/10X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_10X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanomonsv/10X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanomonsv_10X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/nanomonsv/10X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_5X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanomonsv/5X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanomonsv_5X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/minimap2/nanomonsv/5X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"

rule nanomonsv_call_5X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanomonsv/5X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanomonsv_5X_call.log"
    threads: 1
    params:
        wd=f"{RESULTDIR}/ngmlr/nanomonsv/5X",
    shell:
        "mkdir -p {params.wd} && nanomonsv parse {input.bam} {params.wd}/GM24385 && nanomonsv get --min_tumor_variant_read_num 2 --use_racon {params.wd}/GM24385 {input.bam} {input.ref} && mv {params.wd}/GM24385.nanomonsv.result.vcf {output} 2>{log}"
