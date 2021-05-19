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
        for var in ${{supp}}; do bcftools view -f PASS -i "DV>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","CUTESV","minimap2"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","CUTESV","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","CUTESV","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
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
        for var in ${{supp}}; do bcftools view -f PASS -i "DV>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","CUTESV","ngmlr"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","CUTESV","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","CUTESV","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule lra_cutesv_bysupport:
    input:
        testvcf=f"{RESULTDIR}/lra/cutesv/total/GM24385.vcf",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/lra/cutesv/total/cutesv.lra.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/cutesv_lra_prf1bysupport.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/lra/cutesv/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/lra/cutesv/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/lra/cutesv/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/lra/cutesv/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do bcftools view -f PASS -i "DV>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","CUTESV","lra"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","CUTESV","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","CUTESV","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
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
        for var in ${{supp}}; do cat <(cat {input.testvcf}| grep "^#") <(cat {input.testvcf}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | grep -v -E "BND|DUP|INV" | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{var}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","SNIFFLES","minimap2"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","SNIFFLES","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","SNIFFLES","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
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
        for var in ${{supp}}; do cat <(cat {input.testvcf}| grep "^#") <(cat {input.testvcf}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | grep -v -E "BND|DUP|INV" | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{var}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","SNIFFLES","ngmlr"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","SNIFFLES","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","SNIFFLES","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule lra_sniffles_bysupport:
    input:
        testvcf=f"{RESULTDIR}/lra/sniffles/total/GM24385.vcf",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/lra/sniffles/total/sniffles.lra.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/sniffles_lra_prf1bysupport.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/lra/sniffles/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/lra/sniffles/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/lra/sniffles/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/lra/sniffles/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do cat <(cat {input.testvcf}| grep "^#") <(cat {input.testvcf}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | grep -v -E "BND|DUP|INV" | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{var}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","SNIFFLES","lra"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","SNIFFLES","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","SNIFFLES","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
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
        for var in ${{supp}}; do bcftools view -f PASS -i "SUPPORT>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","SVIM","minimap2"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","SVIM","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","SVIM","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
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
        for var in ${{supp}}; do bcftools view -f PASS -i "SUPPORT>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","SVIM","ngmlr"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","SVIM","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","SVIM","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule lra_svim_bysupport:
    input:
        testvcf=f"{RESULTDIR}/lra/svim/total/GM24385.vcf",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/lra/svim/total/svim.lra.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/svim_lra_prf1bysupport.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/lra/svim/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/lra/svim/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/lra/svim/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/lra/svim/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do bcftools view -f PASS -i "SUPPORT>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","SVIM","lra"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","SVIM","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","SVIM","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
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
        for var in ${{supp}}; do bcftools view -f PASS -i "AD[:1]>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} --includebed {input.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","PBSV","minimap2"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","PBSV","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","PBSV","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"PBSV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule combine_results_bysupport:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/{{tool}}/total/{{tool}}.{{aligner}}.bysupport.tsv", tool=["cutesv", "svim", "sniffles"], aligner=["minimap2", "ngmlr", "lra"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/pbsv.pbmm2.bysupport.tsv"
    output:
        f"{RESULTDIR}/GM24385.prf1.bysupport.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/combine_prf1bysupport.log"
    shell:
        """
        awk "FNR>1 || NR==1" {input} > {output} 2>{log}
        """

rule plot_results_bysupport:
    input:
        f"{RESULTDIR}/GM24385.prf1.bysupport.tsv"
    output:
        f"{RESULTDIR}/GM24385.prf1.bysupport.pdf"
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
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "DV>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule ngmlr_cutesv_bycoverage:
    input:
        expand(f"{RESULTDIR}/ngmlr/cutesv/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/ngmlr/cutesv/cutesv.ngmlr.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/cutesv_ngmlr_prf1bycoverage.log"
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
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "DV>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule lra_cutesv_bycoverage:
    input:
        expand(f"{RESULTDIR}/lra/cutesv/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/lra/cutesv/cutesv.lra.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/cutesv_lra_prf1bycoverage.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/lra/cutesv/tmp.vcf",
        vcfgz=f"{RESULTDIR}/lra/cutesv/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/lra/cutesv/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/lra/cutesv",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/GM24385.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "DV>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """


rule minimap2_sniffles_bycoverage:
    input:
        expand(f"{RESULTDIR}/minimap2/sniffles/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/minimap2/sniffles/sniffles.minimap2.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/sniffles_minimap2_prf1bycoverage.log"
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



rule lra_sniffles_bycoverage:
    input:
        expand(f"{RESULTDIR}/lra/sniffles/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/lra/sniffles/sniffles.lra.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/sniffles_lra_prf1bycoverage.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/lra/sniffles/tmp.vcf",
        vcfgz=f"{RESULTDIR}/lra/sniffles/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/lra/sniffles/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/lra/sniffles",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/GM24385.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && cat <(cat ${{var}}| grep "^#") <(cat ${{var}}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | grep -v -E "BND|DUP|INV" | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{s}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule minimap2_svim_bycoverage:
    input:
        expand(f"{RESULTDIR}/minimap2/svim/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/minimap2/svim/svim.minimap2.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/svim_minimap2_prf1bycoverage.log"
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

rule lra_svim_bycoverage:
    input:
        expand(f"{RESULTDIR}/lra/svim/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/lra/svim/svim.lra.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/svim_lra_prf1bycoverage.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/lra/svim/tmp.vcf",
        vcfgz=f"{RESULTDIR}/lra/svim/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/lra/svim/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/lra/svim",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/GM24385.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "SUPPORT>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule pbmm2_pbsv_bycoverage:
    input:
        expand(f"{RESULTDIR}/pbmm2/pbsv/{{coverage}}/GM24385.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/pbmm2/pbsv/pbsv.pbmm2.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/pbsv_pbmm2_prf1bycoverage.log"
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
        expand(f"{RESULTDIR}/{{aligner}}/{{tool}}/{{tool}}.{{aligner}}.bycoverage.tsv", tool=["cutesv", "svim", "sniffles"], aligner=["minimap2", "ngmlr", "lra"]),
        f"{RESULTDIR}/pbmm2/pbsv/pbsv.pbmm2.bycoverage.tsv"
    output:
        f"{RESULTDIR}/GM24385.prf1.bycoverage.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/combine_prf1bycoverage.log"
    shell:
        """
        awk "FNR>1 || NR==1" {input} > {output} 2>{log}
        """

rule plot_results_bycoverage:
    input:
        f"{RESULTDIR}/GM24385.prf1.bycoverage.tsv"
    output:
        f"{RESULTDIR}/GM24385.prf1.bycoverage.pdf"
    threads: 1
    conda: "../envs/plot.yaml"    
    log:
        f"{LOGDIR}/results/plot_prf1bycoverage.log"
    shell:
        "Rscript {SCRIPTDIR}/prf1bycoverage.R {RESULTDIR} 2> {log}"

rule prf1_minimap2_cutesv_sniffles:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv", "sniffles"]),
    output:
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_sniffles.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_sniffles.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpminimap2cutesvsniffles.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/minimap2cutesvsniffles.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_cutesv_sniffles.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_sniffles && cat {params.outdir}/truvari_minimap2_cutesv_sniffles/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_svim:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv", "svim"]),
    output:
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_svim.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_svim.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpminimap2cutesvsvim.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/minimap2cutesvsvim.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_cutesv_svim.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_svim && cat {params.outdir}/truvari_minimap2_cutesv_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-svim","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_pbsv:    
    input:
        f"{RESULTDIR}/minimap2/cutesv/total/GM24385.clean.vcf",
        f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.clean.vcf",
    output:
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_pbsv.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_pbsv.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpminimap2cutesvpbsv.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/minimap2cutesvpbsv.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_cutesv_pbsv.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_pbsv && cat {params.outdir}/truvari_minimap2_cutesv_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_sniffles_pbsv:    
    input:
        f"{RESULTDIR}/minimap2/sniffles/total/GM24385.clean.vcf",
        f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.clean.vcf",
    output:
        f"{RESULTDIR}/calls_combined/minimap2_sniffles_pbsv.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/sniffles_pbsv.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpminimap2snifflespbsv.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/minimap2snifflespbsv.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_sniffles_pbsv.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_sniffles_pbsv && cat {params.outdir}/truvari_minimap2_sniffles_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"sniffles-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_svim_pbsv:
    input:
        f"{RESULTDIR}/minimap2/svim/total/GM24385.clean.vcf",
        f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.clean.vcf",
    output:
        f"{RESULTDIR}/calls_combined/minimap2_svim_pbsv.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/svim_pbsv.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpminimap2svimpbsv.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/minimap2svimpbsv.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_svim_pbsv.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_svim_pbsv && cat {params.outdir}/truvari_minimap2_svim_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"svim-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/GM24385.clean.vcf", tools=["sniffles", "svim"]),
    output:
        f"{RESULTDIR}/calls_combined/minimap2_sniffles_svim.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/sniffles_svim.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpminimap2snifflessvim.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/minimap2snifflessvim.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_sniffles_svim.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_sniffles_svim && cat {params.outdir}/truvari_minimap2_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"sniffles-svim","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv","sniffles", "svim"]),
    output:
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_sniffles_svim.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_sniffles_svim.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpminimap2cutesvsnifflessvim.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/minimap2cutesvsnifflessvim.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_cutesv_sniffles_svim.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_sniffles_svim && cat {params.outdir}/truvari_minimap2_cutesv_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles-svim","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_sniffles_pbsv:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv","sniffles"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.clean.vcf"
    output:
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_sniffles_pbsv.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_sniffles_pbsv.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpminimap2cutesvsnifflespbsv.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/minimap2cutesvsnifflespbsv.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_cutesv_sniffles_pbsv.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_sniffles_pbsv && cat {params.outdir}/truvari_minimap2_cutesv_sniffles_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_svim_pbsv:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv","svim"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.clean.vcf"
    output:
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_svim_pbsv.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_svim_pbsv.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpminimap2cutesvsvimpbsv.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/minimap2cutesvsvimpbsv.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_cutesv_svim_pbsv.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_svim_pbsv && cat {params.outdir}/truvari_minimap2_cutesv_svim_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-svim-pbsv","minimap2"}}' >> {output}
        """ 


rule prf1_minimap2_sniffles_svim_pbsv:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/GM24385.clean.vcf", tools=["sniffles","svim"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.clean.vcf"
    output:
        f"{RESULTDIR}/calls_combined/minimap2_sniffles_svim_pbsv.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/sniffles_svim_pbsv.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpminimap2snifflessvimpbsv.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/minimap2snifflessvimpbsv.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_sniffles_svim_pbsv.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_sniffles_svim_pbsv && cat {params.outdir}/truvari_minimap2_sniffles_svim_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"sniffles-svim-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_sniffles_svim_pbsv:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv","sniffles","svim"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/GM24385.clean.vcf"
    output:
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_sniffles_svim_pbsv.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_sniffles_svim_pbsv.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpminimap2cutesvsnifflessvimpbsv.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/minimap2cutesvsnifflessvimpbsv.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_cutesv_sniffles_svim_pbsv.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 4 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_sniffles_svim_pbsv && cat {params.outdir}/truvari_minimap2_cutesv_sniffles_svim_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles-svim-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_ngmlr_cutesv_sniffles:    
    input:
        expand(f"{RESULTDIR}/ngmlr/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv", "sniffles"]),
    output:
        f"{RESULTDIR}/calls_combined/ngmlr_cutesv_sniffles.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_sniffles.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpngmlrcutesvsniffles.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/ngmlrcutesvsniffles.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/ngmlr_cutesv_sniffles.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_ngmlr_cutesv_sniffles && cat {params.outdir}/truvari_ngmlr_cutesv_sniffles/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles","ngmlr"}}' >> {output}
        """ 

rule prf1_ngmlr_cutesv_svim:    
    input:
        expand(f"{RESULTDIR}/ngmlr/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv", "svim"]),
    output:
        f"{RESULTDIR}/calls_combined/ngmlr_cutesv_svim.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_svim.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpngmlrcutesvsvim.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/ngmlrcutesvsvim.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/ngmlr_cutesv_svim.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_ngmlr_cutesv_svim && cat {params.outdir}/truvari_ngmlr_cutesv_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-svim","ngmlr"}}' >> {output}
        """ 

rule prf1_ngmlr_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/ngmlr/{{tools}}/total/GM24385.clean.vcf", tools=["sniffles", "svim"]),
    output:
        f"{RESULTDIR}/calls_combined/ngmlr_sniffles_svim.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/sniffles_svim.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpngmlrsnifflessvim.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/ngmlrsnifflessvim.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/ngmlr_sniffles_svim.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_ngmlr_sniffles_svim && cat {params.outdir}/truvari_ngmlr_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"sniffles-svim","ngmlr"}}' >> {output}
        """ 

rule prf1_ngmlr_cutesv_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/ngmlr/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv","sniffles", "svim"]),
    output:
        f"{RESULTDIR}/calls_combined/ngmlr_cutesv_sniffles_svim.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_sniffles_svim.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmpngmlrcutesvsnifflessvim.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/ngmlrcutesvsnifflessvim.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/ngmlr_cutesv_sniffles_svim.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_ngmlr_cutesv_sniffles_svim && cat {params.outdir}/truvari_ngmlr_cutesv_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles-svim","ngmlr"}}' >> {output}
        """ 

rule prf1_lra_cutesv_sniffles:    
    input:
        expand(f"{RESULTDIR}/lra/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv", "sniffles"]),
    output:
        f"{RESULTDIR}/calls_combined/lra_cutesv_sniffles.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_sniffles.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmplracutesvsniffles.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/lracutesvsniffles.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/lra_cutesv_sniffles.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_lra_cutesv_sniffles && cat {params.outdir}/truvari_lra_cutesv_sniffles/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles","lra"}}' >> {output}
        """ 

rule prf1_lra_cutesv_svim:    
    input:
        expand(f"{RESULTDIR}/lra/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv", "svim"]),
    output:
        f"{RESULTDIR}/calls_combined/lra_cutesv_svim.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_svim.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmplracutesvsvim.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/lracutesvsvim.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/lra_cutesv_svim.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_lra_cutesv_svim && cat {params.outdir}/truvari_lra_cutesv_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-svim","lra"}}' >> {output}
        """ 

rule prf1_lra_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/lra/{{tools}}/total/GM24385.clean.vcf", tools=["sniffles", "svim"]),
    output:
        f"{RESULTDIR}/calls_combined/lra_sniffles_svim.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/sniffles_svim.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmplrasnifflessvim.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/lrasnifflessvim.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/lra_sniffles_svim.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_lra_sniffles_svim && cat {params.outdir}/truvari_lra_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"sniffles-svim","lra"}}' >> {output}
        """ 

rule prf1_lra_cutesv_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/lra/{{tools}}/total/GM24385.clean.vcf", tools=["cutesv","sniffles", "svim"]),
    output:
        f"{RESULTDIR}/calls_combined/lra_cutesv_sniffles_svim.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/cutesv_sniffles_svim.combined.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/tmplracutesvsnifflessvim.txt",
        mergedvcf=f"{RESULTDIR}/calls_combined/lracutesvsnifflessvim.vcf",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/lra_cutesv_sniffles_svim.vcf.gz",
        truthvcf=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.vcf.gz",
        truthbed=f"{VCFDIR}/HG002_SVs_Tier1_v0.6.bed",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools view -i 'SVTYPE=="DEL" || SVTYPE =="INS"' {params.mergedvcf} | bcftools sort -O z -o {params.mergedvcfgz} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --includebed {params.truthbed} --passonly --giabreport -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_lra_cutesv_sniffles_svim && cat {params.outdir}/truvari_lra_cutesv_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles-svim","lra"}}' >> {output}
        """ 


rule combine_results_bycombo:
    input:
        expand(f"{RESULTDIR}/calls_combined/{{aligner}}_cutesv_sniffles.tsv", aligner=["minimap2","ngmlr","lra"]),
        expand(f"{RESULTDIR}/calls_combined/{{aligner}}_cutesv_svim.tsv", aligner=["minimap2","ngmlr","lra"]),
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_pbsv.tsv",
        f"{RESULTDIR}/calls_combined/minimap2_sniffles_pbsv.tsv",
        f"{RESULTDIR}/calls_combined/minimap2_svim_pbsv.tsv",
        expand(f"{RESULTDIR}/calls_combined/{{aligner}}_sniffles_svim.tsv",aligner=["minimap2","ngmlr","lra"]),
        expand(f"{RESULTDIR}/calls_combined/{{aligner}}_cutesv_sniffles_svim.tsv",aligner=["minimap2","ngmlr","lra"]),
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_sniffles_pbsv.tsv",
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_svim_pbsv.tsv",
        f"{RESULTDIR}/calls_combined/minimap2_sniffles_svim_pbsv.tsv",
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_sniffles_svim_pbsv.tsv"
    output:
        f"{RESULTDIR}/GM24385.prf1.bycombo.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/combine_prf1bycombo.log"
    shell:
        """
        awk "FNR>1 || NR==1" {input} > {output} 2>{log}
        """

rule plot_results_bycombo:
    input:
        f"{RESULTDIR}/GM24385.prf1.bycombo.tsv"
    output:
        f"{RESULTDIR}/GM24385.prf1.bycombo.pdf"
    threads: 1
    conda: "../envs/plot.yaml"    
    log:
        f"{LOGDIR}/results/plot_prf1bycombo.log"
    shell:
        "Rscript {SCRIPTDIR}/prf1bycombo.R {RESULTDIR} 2> {log}"


rule combine_stats_by_length:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/{{tool}}/total/{{tool}}.{{aligner}}.bysupport.tsv", tool=["cutesv", "svim", "sniffles"], aligner=["minimap2", "ngmlr","lra"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/pbsv.pbmm2.bysupport.tsv"
    output:
        f"{RESULTDIR}/GM24385.tpfpfn.bylength.tsv"
    log:
        f"{LOGDIR}/results/combine_stats_by_length.log"
    shell:
        """
        echo -e "size\tsvtype\tvartype\ttool\taligner" > {output} && \
        find {RESULTDIR} -wholename '*truvari_supporting_10/svlen.info.tsv*' -exec cat {{}} + >> {output} 2> {log}
        """


rule plot_stats_by_length:
    input:
        f"{RESULTDIR}/GM24385.tpfpfn.bylength.tsv"
    output:
        f"{RESULTDIR}/GM24385.tpfpfn.bylength.pdf",
        f"{RESULTDIR}/GM24385.svtype.bylength.pdf",
        f"{RESULTDIR}/GM24385.prf1.bysvtype.pdf"
    threads: 1
    conda: "../envs/plot.yaml"    
    log:
        f"{LOGDIR}/results/plot_stats_by_length"
    shell:
        "Rscript {SCRIPTDIR}/tpfpfnbysize.R {RESULTDIR} 2> {log}"




