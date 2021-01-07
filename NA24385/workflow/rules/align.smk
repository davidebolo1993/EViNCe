rule minimap2_align:
    input:
       config["genome"],
       config["samples"].values()
    output:
        f"{ALIGNDIR}/GM24385.minimap2.srt.bam"
    threads: 20
    log:
        f"{LOGDIR}/alignments/minimap2_samtools_alignment.log"
    conda: "../envs/minimap2.yaml"
    params:
        rg=r"'@RG\tID:nanopore\tSM:GM24385'"
    shell:
        "minimap2 --MD -ax map-ont -t {threads} -R {params.rg} {input} | samtools sort -@ {threads} -o {output} - 2> {log}"

rule ngmlr_align:
    input:
       ref=config["genome"],
       fq=config["samples"].values()
    output:
        f"{ALIGNDIR}/GM24385.ngmlr.srt.bam"
    threads: 20
    log:
        f"{LOGDIR}/alignments/ngmlr_samtools_alignment.log"
    conda: "../envs/ngmlr.yaml"
    shell:
        "zcat {input.fq} | ngmlr -r {input.ref} -x ont -t {threads} --rg-id nanopore --rg-sm GM24385  | samtools sort -@ {threads} -o {output} - 2> {log}"

rule samtools_index_minimap2:
    input:
        f"{ALIGNDIR}/GM24385.minimap2.srt.bam"
    output:
        f"{ALIGNDIR}/GM24385.minimap2.srt.bam.bai"
    conda: "../envs/minimap2.yaml"
    threads: 10
    shell:
        "samtools index -@ {threads} {input}"

rule samtools_index_ngmlr:
    input:
        f"{ALIGNDIR}/GM24385.ngmlr.srt.bam"
    output:
        f"{ALIGNDIR}/GM24385.ngmlr.srt.bam.bai"
    conda: "../envs/minimap2.yaml"
    threads: 10
    shell:
        "samtools index -@ {threads} {input}"
