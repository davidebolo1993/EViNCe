rule minimap2_cutesv_bysupport:
    input:
        testvcf=f"{RESULTDIR}/minimap2/cutesv/total/GM24385.vcf",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/minimap2/cutesv/total/cutesv.minimap2.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/cutesv_minimap2_prf1bysupport.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/minimap2/cutesv/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/minimap2/cutesv/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/minimap2/cutesv/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/minimap2/cutesv/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do bcftools view -f PASS -i "DV>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -e 'SVTYPE=="BND" || SVTYPE =="DUP" || SVTYPE=="INV"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule ngmlr_cutesv_bysupport:
    input:
        testvcf=f"{RESULTDIR}/ngmlr/cutesv/total/GM24385.vcf",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/ngmlr/cutesv/total/cutesv.ngmlr.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/cutesv_ngmlr_prf1bysupport.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/ngmlr/cutesv/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/ngmlr/cutesv/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/ngmlr/cutesv/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/ngmlr/cutesv/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do bcftools view -f PASS -i "DV>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -e 'SVTYPE=="BND" || SVTYPE =="DUP" || SVTYPE=="INV"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule minimap2_sniffles_bysupport:
    input:
        testvcf=f"{RESULTDIR}/minimap2/sniffles/total/GM24385.vcf",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/minimap2/sniffles/total/sniffles.minimap2.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/sniffles_minimap2_prf1bysupport.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/minimap2/sniffles/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/minimap2/sniffles/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/minimap2/sniffles/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/minimap2/sniffles/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do cat <(cat {input.testvcf}| grep "^#") <(cat {input.testvcf}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | grep -v -E "BND|DUP|INV" | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{var}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule ngmlr_sniffles_bysupport:
    input:
        testvcf=f"{RESULTDIR}/ngmlr/sniffles/total/GM24385.vcf",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/ngmlr/sniffles/total/sniffles.ngmlr.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/sniffles_ngmlr_prf1bysupport.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/ngmlr/sniffles/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/ngmlr/sniffles/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/ngmlr/sniffles/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/ngmlr/sniffles/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do cat <(cat {input.testvcf}| grep "^#") <(cat {input.testvcf}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | grep -v -E "BND|DUP|INV" | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{var}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule minimap2_svim_bysupport:
    input:
        testvcf=f"{RESULTDIR}/minimap2/svim/total/GM24385.vcf",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/minimap2/svim/total/svim.minimap2.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/svim_minimap2_prf1bysupport.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/minimap2/svim/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/minimap2/svim/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/minimap2/svim/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/minimap2/svim/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do bcftools view -f PASS -i "SUPPORT>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule ngmlr_svim_bysupport:
    input:
        testvcf=f"{RESULTDIR}/ngmlr/svim/total/GM24385.vcf",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/ngmlr/svim/total/svim.ngmlr.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/svim_ngmlr_prf1bysupport.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/ngmlr/svim/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/ngmlr/svim/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/ngmlr/svim/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/ngmlr/svim/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do bcftools view -f PASS -i "SUPPORT>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule pbmm2_pbsv_bysupport:
    input:
        testvcf=f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.vcf",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/pbmm2/pbsv/total/pbsv.pbmm2.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/svim_minimap2_prf1bysupport.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/pbmm2/pbsv/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/pbmm2/pbsv/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/pbmm2/pbsv/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/pbmm2/pbsv/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do bcftools view -f PASS -i "AD[:1]>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"PBSV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule combine_results_bysupport:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/{{tool}}/total/{{tool}}.{{aligner}}.bysupport.tsv", tool=["cutesv", "svim", "sniffles"], aligner=["minimap2", "ngmlr"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/pbsv.pbmm2.bysupport.tsv"
    output:
        f"{RESULTDIR}/prf1.bysupport.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/combine_prf1bysupport.log"
    shell:
        """
        awk "FNR>1 || NR==1" {input} > {output} 2>{log}
        """

rule plot_results_bysupport:
    input:
        f"{RESULTDIR}/prf1.bysupport.tsv"
    output:
        f"{RESULTDIR}/prf1.bysupport.pdf"
    threads: 1
    conda: "../envs/plot.yaml"    
    log:
        f"{LOGDIR}/results/plot_prf1bysupport.log"
    shell:
        "Rscript {SCRIPTDIR}/prf1bysupport.R {RESULTDIR} 2> {log}"


rule minimap2_cutesv_bycoverage:
    input:
        expand(f"{RESULTDIR}/minimap2/cutesv/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/minimap2/cutesv/cutesv.minimap2.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/cutesv_minimap2_prf1bycoverage.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/minimap2/cutesv/tmp.vcf",
        vcfgz=f"{RESULTDIR}/minimap2/cutesv/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/minimap2/cutesv/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/minimap2/cutesv",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/GM24385.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "DV>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -e 'SVTYPE=="BND" || SVTYPE =="DUP" || SVTYPE=="INV"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """


rule ngmlr_cutesv_bycoverage:
    input:
        expand(f"{RESULTDIR}/ngmlr/cutesv/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/ngmlr/cutesv/cutesv.ngmlr.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/cutesv_ngmlr_prf1bycoverage.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/ngmlr/cutesv/tmp.vcf",
        vcfgz=f"{RESULTDIR}/ngmlr/cutesv/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/ngmlr/cutesv/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/ngmlr/cutesv",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/GM24385.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "DV>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -e 'SVTYPE=="BND" || SVTYPE =="DUP" || SVTYPE=="INV"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule minimap2_sniffles_bycoverage:
    input:
        expand(f"{RESULTDIR}/minimap2/sniffles/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/minimap2/sniffles/sniffles.minimap2.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/sniffles_minimap2_prf1bycoverage.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/minimap2/sniffles/tmp.vcf",
        vcfgz=f"{RESULTDIR}/minimap2/sniffles/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/minimap2/sniffles/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/minimap2/sniffles",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/GM24385.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && cat <(cat ${{var}}| grep "^#") <(cat ${{var}}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | grep -v -E "BND|DUP|INV" | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{s}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """


rule ngmlr_sniffles_bycoverage:
    input:
        expand(f"{RESULTDIR}/ngmlr/sniffles/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/ngmlr/sniffles/sniffles.ngmlr.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/sniffles_ngmlr_prf1bycoverage.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/ngmlr/sniffles/tmp.vcf",
        vcfgz=f"{RESULTDIR}/ngmlr/sniffles/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/ngmlr/sniffles/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/ngmlr/sniffles",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/GM24385.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && cat <(cat ${{var}}| grep "^#") <(cat ${{var}}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | grep -v -E "BND|DUP|INV" | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{s}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """


rule minimap2_svim_bycoverage:
    input:
        expand(f"{RESULTDIR}/minimap2/svim/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/minimap2/svim/svim.minimap2.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/svim_minimap2_prf1bycoverage.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/minimap2/svim/tmp.vcf",
        vcfgz=f"{RESULTDIR}/minimap2/svim/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/minimap2/svim/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/minimap2/svim",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/GM24385.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "SUPPORT>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule ngmlr_svim_bycoverage:
    input:
        expand(f"{RESULTDIR}/ngmlr/svim/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/ngmlr/svim/svim.ngmlr.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/svim_ngmlr_prf1bycoverage.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/ngmlr/svim/tmp.vcf",
        vcfgz=f"{RESULTDIR}/ngmlr/svim/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/ngmlr/svim/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/ngmlr/svim",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/GM24385.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "SUPPORT>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule pbmm2_pbsv_bycoverage:
    input:
        expand(f"{RESULTDIR}/pbmm2/pbsv/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/pbmm2/pbsv/pbsv.pbmm2.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/pbsv_pbmm2_prf1bycoverage.log"
    conda: "../envs/sumsvs.yaml" #truvari must be installed manually as this is outdated in bioconda
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/pbmm2/pbsv/tmp.vcf",
        vcfgz=f"{RESULTDIR}/pbmm2/pbsv/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/pbmm2/pbsv/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/pbmm2/pbsv",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/GM24385.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "AD[:1]>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"PBSV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule combine_results_bycoverage:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/{{tool}}/{{tool}}.{{aligner}}.bycoverage.tsv", tool=["cutesv", "svim", "sniffles"], aligner=["minimap2", "ngmlr"]),
        f"{RESULTDIR}/pbmm2/pbsv/pbsv.pbmm2.bycoverage.tsv"
    output:
        f"{RESULTDIR}/prf1.bycoverage.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/combine_prf1bycoverage.log"
    shell:
        """
        awk "FNR>1 || NR==1" {input} > {output} 2>{log}
        """

rule plot_results_bycoverage:
    input:
        f"{RESULTDIR}/prf1.bycoverage.tsv"
    output:
        f"{RESULTDIR}/prf1.bycoverage.pdf"
    threads: 1
    conda: "../envs/plot.yaml"    
    log:
        f"{LOGDIR}/results/plot_prf1bycoverage.log"
    shell:
        "Rscript {SCRIPTDIR}/prf1bycoverage.R {RESULTDIR} 2> {log}"


