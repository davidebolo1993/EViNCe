rule cutesv_call_total_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/total/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_total_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/total/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_total_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/total/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_total_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/total/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_5X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/5X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_5X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/5X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_5X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/5X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/5X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_5X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/5X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_10X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/10X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_10X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/10X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_10X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/10X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/10X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_10X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/10X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}"


rule cutesv_call_15X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/15X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_15X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/15X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_15X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/15X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/15X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_15X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/15X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_20X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/20X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_20X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/20X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_20X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/20X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/20X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_20X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/20X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}"


rule cutesv_call_25X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/25X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_25X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/25X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_25X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/25X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/25X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_25X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/25X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}"


rule cutesv_call_35X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/GM24385.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/35X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_35X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/35X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_35X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/GM24385.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/35X/GM24385.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/35X/GM24385.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_35X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/35X/tmpdir"
    shell:
        "cuteSV -s 3 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -S GM24385 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

