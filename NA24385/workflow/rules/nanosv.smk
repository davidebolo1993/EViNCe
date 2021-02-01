rule nanosv_call_total_minimap2:
    input:
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanosv/total/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanosv_total_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/nanosv/total",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}" 

rule nanosv_call_total_ngmlr:
    input:
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanosv/total/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanosv_total_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/nanosv/total",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}"

rule nanosv_call_5X_minimap2:
    input:
        bam=f"{ALIGNDIR}/5X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanosv/5X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanosv_5X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/nanosv/5X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}" 

rule nanosv_call_5X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/5X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanosv/5X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanosv_5X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/nanosv/5X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}"

rule nanosv_call_10X_minimap2:
    input:
        bam=f"{ALIGNDIR}/10X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanosv/10X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanosv_10X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/nanosv/10X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}" 

rule nanosv_call_10X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/10X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanosv/10X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanosv_10X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/nanosv/10X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}"

rule nanosv_call_15X_minimap2:
    input:
        bam=f"{ALIGNDIR}/15X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanosv/15X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanosv_15X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/nanosv/15X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}" 

rule nanosv_call_15X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/15X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanosv/15X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanosv_15X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/nanosv/15X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}"

rule nanosv_call_20X_minimap2:
    input:
        bam=f"{ALIGNDIR}/20X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanosv/20X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanosv_20X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/nanosv/20X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}" 

rule nanosv_call_20X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/20X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanosv/20X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanosv_20X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/nanosv/20X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}"


rule nanosv_call_25X_minimap2:
    input:
        bam=f"{ALIGNDIR}/25X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanosv/25X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanosv_25X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/nanosv/25X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}" 

rule nanosv_call_25X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/25X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanosv/25X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanosv_25X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/nanosv/25X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}"


rule nanosv_call_35X_minimap2:
    input:
        bam=f"{ALIGNDIR}/35X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/nanosv/35X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_nanosv_35X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/nanosv/35X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}" 

rule nanosv_call_35X_ngmlr:
    input:
        bam=f"{ALIGNDIR}/35X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/nanosv/35X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_nanosv_35X_call.log"
    conda: "../envs/nanosv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/nanosv/35X",
        bed=config["nanosvbed"]
    shell:
        "mkdir -p {params.wd} && NanoSV -t {threads} -o {output} -b {params.bed} -s samtools {input.bam} 2>{log}"
