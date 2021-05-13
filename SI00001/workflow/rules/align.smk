rule minimap2_align:
    input:
       config["genome"],
       config["samples"].values()
    output:
        f"{ALIGNDIR}/SI00001.minimap2.srt.bam"
    threads: 20
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_alignment.log"
    conda: "../envs/minimap2.yaml"
    params:
        rg=r"'@RG\tID:nanopore\tSM:SI00001'"
    shell:
        "minimap2 --MD -ax map-ont -t {threads} -R {params.rg} {input} | samtools sort -@ {threads} -o {output} - 2>{log}"

rule ngmlr_align:
    input:
       ref=config["genome"],
       fq=config["samples"].values()
    output:
        f"{ALIGNDIR}/SI00001.ngmlr.srt.bam"
    threads: 20
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_alignment.log"
    conda: "../envs/ngmlr.yaml"
    shell:
        "cat {input.fq} | ngmlr -r {input.ref} -x ont -t {threads} --rg-id nanopore --rg-sm SI00001  | samtools sort -@ {threads} -o {output} - 2>{log}"

rule lra_align:
    input:
       ref=config["genome"],
       fq=config["samples"].values()
    output:
        f"{ALIGNDIR}/SI00001.lra.srt.bam"
    threads: 20
    log:
        f"{LOGDIR}/alignments/lra_samtools_alignment.log"
    conda: "../envs/lra.yaml"
    params:
        rg=r"'@RG\tID:nanopore\tSM:SI00001'",
        tmpbam=f"{ALIGNDIR}/SI00001.lra.srt.tmp.bam"
    shell:
        "lra index -ONT {input.ref} && cat {input.fq} | lra align -ONT {input.ref} /dev/stdin -t {threads} -p s | samtools addreplacerg -@ {threads} -r {params.rg} - | samtools sort -@ {threads} -o {params.tmpbam} - && samtools calmd -@ {threads} -b {params.tmpbam} {input.ref} > {output} && rm {params.tmpbam} 2>{log}"

rule pbmm2_align_index:
    input:
       config["genome"],
       config["fofn"]
    output:
        bam=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam",
        bai=f"{ALIGNDIR}/SI00001.pbmm2.srt.bam.bai"
    threads: 20
    log:
        f"{LOGDIR}/alignments/pbmm2_alignment.log"
    conda: "../envs/pbmm2.yaml"
    params:
        rg=r"'@RG\tID:nanopore\tSM:SI00001'"
    shell:
        "pbmm2 align --preset subread --sort --bam-index BAI --rg {params.rg} -j {threads} {input} {output.bam} 2>{log}"
        
rule samtools_index_minimap2:
    input:
        f"{ALIGNDIR}/SI00001.minimap2.srt.bam"
    output:
        f"{ALIGNDIR}/SI00001.minimap2.srt.bam.bai"
    conda: "../envs/minimap2.yaml"
    threads: 10
    shell:
        "samtools index -@ {threads} {input}"

rule samtools_index_lra:
    input:
        f"{ALIGNDIR}/SI00001.lra.srt.bam"
    output:
        f"{ALIGNDIR}/SI00001.lra.srt.bam.bai"
    conda: "../envs/lra.yaml"
    threads: 10
    shell:
        "samtools index -@ {threads} {input}"

rule samtools_index_ngmlr:
    input:
        f"{ALIGNDIR}/SI00001.ngmlr.srt.bam"
    output:
        f"{ALIGNDIR}/SI00001.ngmlr.srt.bam.bai"
    conda: "../envs/ngmlr.yaml"
    threads: 10
    shell:
        "samtools index -@ {threads} {input}"



