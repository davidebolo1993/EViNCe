rule pbsv_call_total_pbmm2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.pbmm2.srt.bam.bai"
    output:
        f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.vcf"
    log:
        f"{LOGDIR}/results/pbmm2_pbsv_total_call.log"
    conda: "../envs/pbsv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/pbmm2/pbsv/total",
        sigs=f"{RESULTDIR}/pbmm2/pbsv/total/signatures.svsig.gz"
    shell:
        "mkdir -p {params.wd} && pbsv discover {input.bam} {params.sigs} 2>{log} &&  pbsv call -j {threads} -m 50 {input.ref} {params.sigs} {output} 2>>{log}" 

rule pbsv_call_5X_pbmm2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.pbmm2.srt.bam.bai"
    output:
        f"{RESULTDIR}/pbmm2/pbsv/5X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/pbmm2_pbsv_5X_call.log"
    conda: "../envs/pbsv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/pbmm2/pbsv/5X",
        sigs=f"{RESULTDIR}/pbmm2/pbsv/5X/signatures.svsig.gz"
    shell:
        "mkdir -p {params.wd} && pbsv discover {input.bam} {params.sigs} 2>{log} &&  pbsv call -j {threads} -m 50 {input.ref} {params.sigs} {output} 2>>{log}" 

rule pbsv_call_10X_pbmm2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.pbmm2.srt.bam.bai"
    output:
        f"{RESULTDIR}/pbmm2/pbsv/10X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/pbmm2_pbsv_10X_call.log"
    conda: "../envs/pbsv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/pbmm2/pbsv/10X",
        sigs=f"{RESULTDIR}/pbmm2/pbsv/10X/signatures.svsig.gz"
    shell:
        "mkdir -p {params.wd} && pbsv discover {input.bam} {params.sigs} 2>{log} &&  pbsv call -j {threads} -m 50 {input.ref} {params.sigs} {output} 2>>{log}" 


rule pbsv_call_15X_pbmm2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.pbmm2.srt.bam.bai"
    output:
        f"{RESULTDIR}/pbmm2/pbsv/15X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/pbmm2_pbsv_15X_call.log"
    conda: "../envs/pbsv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/pbmm2/pbsv/15X",
        sigs=f"{RESULTDIR}/pbmm2/pbsv/15X/signatures.svsig.gz"
    shell:
        "mkdir -p {params.wd} && pbsv discover {input.bam} {params.sigs} 2>{log} &&  pbsv call -j {threads} -m 50 {input.ref} {params.sigs} {output} 2>>{log}" 


rule pbsv_call_20X_pbmm2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.pbmm2.srt.bam.bai"
    output:
        f"{RESULTDIR}/pbmm2/pbsv/20X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/pbmm2_pbsv_20X_call.log"
    conda: "../envs/pbsv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/pbmm2/pbsv/20X",
        sigs=f"{RESULTDIR}/pbmm2/pbsv/20X/signatures.svsig.gz"
    shell:
        "mkdir -p {params.wd} && pbsv discover {input.bam} {params.sigs} 2>{log} &&  pbsv call -j {threads} -m 50 {input.ref} {params.sigs} {output} 2>>{log}" 


rule pbsv_call_25X_pbmm2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.pbmm2.srt.bam.bai"
    output:
        f"{RESULTDIR}/pbmm2/pbsv/25X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/pbmm2_pbsv_25X_call.log"
    conda: "../envs/pbsv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/pbmm2/pbsv/25X",
        sigs=f"{RESULTDIR}/pbmm2/pbsv/25X/signatures.svsig.gz"
    shell:
        "mkdir -p {params.wd} && pbsv discover {input.bam} {params.sigs} 2>{log} &&  pbsv call -j {threads} -m 50 {input.ref} {params.sigs} {output} 2>>{log}" 


rule pbsv_call_35X_pbmm2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/GM24385.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.pbmm2.srt.bam.bai"
    output:
        f"{RESULTDIR}/pbmm2/pbsv/35X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/pbmm2_pbsv_35X_call.log"
    conda: "../envs/pbsv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/pbmm2/pbsv/35X",
        sigs=f"{RESULTDIR}/pbmm2/pbsv/35X/signatures.svsig.gz"
    shell:
        "mkdir -p {params.wd} && pbsv discover {input.bam} {params.sigs} 2>{log} &&  pbsv call -j {threads} -m 50 {input.ref} {params.sigs} {output} 2>>{log}" 


