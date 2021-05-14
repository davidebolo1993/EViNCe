rule cutesv_call_total_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_total_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/total/tmpdir",
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_total_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_total_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/total/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_total_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/cutesv/total/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_cutesv_total_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/lra/cutesv/total/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_5X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_5X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/5X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_5X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_5X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/5X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_5X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/5X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/5X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/cutesv/5X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_cutesv_5X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/lra/cutesv/5X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_10X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_10X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/10X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_10X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_10X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/10X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_10X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/10X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/10X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/cutesv/10X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_cutesv_10X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/lra/cutesv/10X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_15X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_15X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/15X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_15X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_15X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/15X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"


rule cutesv_call_15X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/15X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/15X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/cutesv/15X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_cutesv_15X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/lra/cutesv/15X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_20X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_20X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/20X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_20X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_20X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/20X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_20X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/20X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/20X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/cutesv/20X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_cutesv_20X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/lra/cutesv/20X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"


rule cutesv_call_25X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_25X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/25X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_25X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_25X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/25X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"


rule cutesv_call_25X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/25X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/25X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/cutesv/25X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_cutesv_25X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/lra/cutesv/25X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_35X_minimap2:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/SI00001.minimap2.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.minimap2.srt.bam.bai"
    output:
        f"{RESULTDIR}/minimap2/cutesv/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/minimap2_cutesv_35X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/minimap2/cutesv/35X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}" 

rule cutesv_call_35X_ngmlr:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/SI00001.ngmlr.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.ngmlr.srt.bam.bai"
    output:
        f"{RESULTDIR}/ngmlr/cutesv/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/ngmlr_cutesv_35X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/ngmlr/cutesv/35X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"

rule cutesv_call_35X_lra:
    input:
        ref=config["genome"],
        bam=f"{ALIGNDIR}/35X/SI00001.lra.srt.bam",
        bai=f"{ALIGNDIR}/35X/SI00001.lra.srt.bam.bai"
    output:
        f"{RESULTDIR}/lra/cutesv/35X/SI00001.vcf"
    log:
        f"{LOGDIR}/results/lra_cutesv_35X_call.log"
    conda: "../envs/cutesv.yaml"
    threads: 10
    params:
        wd=f"{RESULTDIR}/lra/cutesv/35X/tmpdir"
    shell:
        "mkdir -p {params.wd} && cuteSV -s 2 -l 50 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype -t {threads} --retain_work_dir -S SI00001 {input.bam} {input.ref} {output} {params.wd} 2>{log}"



