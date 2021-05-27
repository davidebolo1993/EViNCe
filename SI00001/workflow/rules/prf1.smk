rule minimap2_cutesv_bysupport:
    input:
        testvcf=f"{RESULTDIR}/minimap2/cutesv/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
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
        for var in ${{supp}}; do bcftools view -f PASS -i "DV>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","CUTESV","minimap2"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","CUTESV","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","CUTESV","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule ngmlr_cutesv_bysupport:
    input:
        testvcf=f"{RESULTDIR}/ngmlr/cutesv/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
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
        for var in ${{supp}}; do bcftools view -f PASS -i "DV>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","CUTESV","ngmlr"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","CUTESV","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","CUTESV","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule lra_cutesv_bysupport:
    input:
        testvcf=f"{RESULTDIR}/lra/cutesv/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
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
        for var in ${{supp}}; do bcftools view -f PASS -i "DV>=${{var}}" {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","CUTESV","lra"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","CUTESV","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","CUTESV","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule minimap2_sniffles_bysupport:
    input:
        testvcf=f"{RESULTDIR}/minimap2/sniffles/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
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
        for var in ${{supp}}; do cat <(cat {input.testvcf}| grep "^#") <(cat {input.testvcf}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{var}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","SNIFFLES","minimap2"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","SNIFFLES","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","SNIFFLES","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule ngmlr_sniffles_bysupport:
    input:
        testvcf=f"{RESULTDIR}/ngmlr/sniffles/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
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
        for var in ${{supp}}; do cat <(cat {input.testvcf}| grep "^#") <(cat {input.testvcf}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{var}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","SNIFFLES","ngmlr"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","SNIFFLES","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","SNIFFLES","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule lra_sniffles_bysupport:
    input:
        testvcf=f"{RESULTDIR}/lra/sniffles/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
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
        for var in ${{supp}}; do cat <(cat {input.testvcf}| grep "^#") <(cat {input.testvcf}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{var}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"TP","SNIFFLES","lra"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FP","SNIFFLES","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","SNIFFLES","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule minimap2_svim_bysupport:
    input:
        testvcf=f"{RESULTDIR}/minimap2/svim/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
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
        for var in ${{supp}}; do bcftools view -f PASS -i "SUPPORT>=${{var}}"  {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} |  sed 's/DUP:TANDEM/DUP/g' | sed 's/DUP:INT/DUP/g' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","TP","SVIM","minimap2"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"TP", "SVIM", "minimap2"; else if ($4 == "INV") print $2-$1,"INV","TP", "SVIM", "minimap2"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","FP","SVIM","minimap2"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"FP", "SVIM", "minimap2"; else if ($4 == "INV") print $2-$1,"INV","FP", "SVIM", "minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","FN","SVIM","minimap2"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"FN", "SVIM", "minimap2"; else if ($4 == "INV") print $2-$1,"INV","FN", "SVIM", "minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule ngmlr_svim_bysupport:
    input:
        testvcf=f"{RESULTDIR}/ngmlr/svim/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
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
        for var in ${{supp}}; do bcftools view -f PASS -i "SUPPORT>=${{var}}"  {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} |  sed 's/DUP:TANDEM/DUP/g' | sed 's/DUP:INT/DUP/g' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","TP","SVIM","ngmlr"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"TP", "SVIM", "ngmlr"; else if ($4 == "INV") print $2-$1,"INV","TP", "SVIM", "ngmlr"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","FP","SVIM","ngmlr"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"FP", "SVIM", "ngmlr"; else if ($4 == "INV") print $2-$1,"INV","FP", "SVIM", "ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","FN","SVIM","ngmlr"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"FN", "SVIM", "ngmlr"; else if ($4 == "INV") print $2-$1,"INV","FN", "SVIM", "ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule lra_svim_bysupport:
    input:
        testvcf=f"{RESULTDIR}/lra/svim/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
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
        for var in ${{supp}}; do bcftools view -f PASS -i "SUPPORT>=${{var}}"  {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} |  sed 's/DUP:TANDEM/DUP/g' | sed 's/DUP:INT/DUP/g' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","TP","SVIM","lra"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"TP", "SVIM", "lra"; else if ($4 == "INV") print $2-$1,"INV","TP", "SVIM", "lra"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","FP","SVIM","lra"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"FP", "SVIM", "lra"; else if ($4 == "INV") print $2-$1,"INV","FP", "SVIM", "lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","FN","SVIM","lra"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"FN", "SVIM", "lra"; else if ($4 == "INV") print $2-$1,"INV","FN", "SVIM", "lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """


rule pbmm2_pbsv_bysupport:
    input:
        testvcf=f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/pbmm2/pbsv/total/pbsv.pbmm2.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/pbmm2_pbsv_prf1bysupport.log"
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
        for var in ${{supp}}; do bcftools view -f PASS -i "AD[:1]>=${{var}}"  {input.testvcf} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","TP","PBSV","minimap2"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"TP", "PBSV", "minimap2"; else if ($4 == "INV") print $2-$1,"INV","TP", "PBSV", "minimap2"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","FP","PBSV","minimap2"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"FP", "PBSV", "minimap2"; else if ($4 == "INV") print $2-$1,"INV","FP", "PBSV", "minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv &&  bcftools query -f '%POS\t%INFO/END\t%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{if ($4 == "DEL") print -$3,"DEL","FN","PBSV","minimap2"; else if ($4 == "INS" || $4 == "DUP") print $3,$4,"FN", "PBSV", "minimap2"; else if ($4 == "INV") print $2-$1,"INV","FN", "PBSV", "minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"PBSV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule minimap2_npinv_bysupport:
    input:
        testvcf=f"{RESULTDIR}/minimap2/npinv/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/minimap2/npinv/total/npinv.minimap2.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/npinv_minimap2_prf1bysupport.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/minimap2/npinv/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/minimap2/npinv/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/minimap2/npinv/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/minimap2/npinv/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do cat <(grep "^#" {input.testvcf}) <(grep -v -e "#" -e "_K" -e "HLA" -e "chrUn_" -e "_GL" -e "JH" -e "EB" -e "chrM" {input.testvcf}) | bcftools view -f PASS -i "DV>=${{var}}" | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%POS\t%INFO/END\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print $2-$1,"INV","TP","NPINV","minimap2"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%POS\t%INFO/END\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print $2-$1,"INV","FP","NPINV","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","NPINV","minimap2"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"NPINV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule ngmlr_npinv_bysupport:
    input:
        testvcf=f"{RESULTDIR}/ngmlr/npinv/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/ngmlr/npinv/total/npinv.ngmlr.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/npinv_ngmlr_prf1bysupport.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/ngmlr/npinv/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/ngmlr/npinv/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/ngmlr/npinv/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/ngmlr/npinv/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do cat <(grep "^#" {input.testvcf}) <(grep -v -e "#" -e "_K" -e "HLA" -e "chrUn_" -e "_GL" -e "JH" -e "EB" -e "chrM" {input.testvcf}) | bcftools view -f PASS -i "DV>=${{var}}" | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%POS\t%INFO/END\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print $2-$1,"INV","TP","NPINV","ngmlr"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%POS\t%INFO/END\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print $2-$1,"INV","FP","NPINV","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","NPINV","ngmlr"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"NPINV","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule lra_npinv_bysupport:
    input:
        testvcf=f"{RESULTDIR}/lra/npinv/total/SI00001.vcf",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    output:
        f"{RESULTDIR}/lra/npinv/total/npinv.lra.bysupport.tsv"  
    log:
        f"{LOGDIR}/results/npinv_lra_prf1bysupport.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/lra/npinv/total/tmp.vcf",
        vcfgz=f"{RESULTDIR}/lra/npinv/total/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/lra/npinv/total/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/lra/npinv/total"
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        supp="2 5 10 15 20 25 30 35 40 45 50" && \
        for var in ${{supp}}; do cat <(grep "^#" {input.testvcf}) <(grep -v -e "#" -e "_K" -e "HLA" -e "chrUn_" -e "_GL" -e "JH" -e "EB" -e "chrM" {input.testvcf}) | bcftools view -f PASS -i "DV>=${{var}}" | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {input.ref} -b {input.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_supporting_${{var}} && bcftools query -f '%POS\t%INFO/END\n' {params.outdir}/truvari_supporting_${{var}}/tp-call.vcf | awk '{{OFS=FS="\t"}}{{print $2-$1,"INV","TP","NPINV","lra"}}' > {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%POS\t%INFO/END\n' {params.outdir}/truvari_supporting_${{var}}/fp.vcf | awk '{{OFS=FS="\t"}}{{print $2-$1,"INV","FP","NPINV","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && bcftools query -f '%SVLEN\t%SVTYPE\n' {params.outdir}/truvari_supporting_${{var}}/fn.vcf | awk '{{OFS=FS="\t"}}{{print ($1>0)?$1:-$1,$2,"FN","NPINV","lra"}}' >> {params.outdir}/truvari_supporting_${{var}}/svlen.info.tsv && cat {params.outdir}/truvari_supporting_${{var}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{var}} '{{OFS=FS="\t"}}{{print $0,su,"NPINV","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log}; done
        """

rule combine_results_bysupport:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/{{tool}}/total/{{tool}}.{{aligner}}.bysupport.tsv", tool=["cutesv", "svim", "sniffles", "npinv"], aligner=["minimap2", "ngmlr", "lra"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/pbsv.pbmm2.bysupport.tsv"
    output:
        f"{RESULTDIR}/SI00001.prf1.bysupport.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/combine_prf1bysupport.log"
    shell:
        """
        awk "FNR>1 || NR==1" {input} > {output} 2>{log}
        """

rule plot_results_bysupport:
    input:
        f"{RESULTDIR}/SI00001.prf1.bysupport.tsv"
    output:
        f"{RESULTDIR}/SI00001.prf1.bysupport.pdf"
    threads: 1
    conda: "../envs/plot.yaml"    
    log:
        f"{LOGDIR}/results/plot_prf1bysupport.log"
    shell:
        "Rscript {SCRIPTDIR}/prf1bysupport.R {RESULTDIR} 2> {log}"

rule minimap2_cutesv_bycoverage:
    input:
        expand(f"{RESULTDIR}/minimap2/cutesv/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "DV>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule ngmlr_cutesv_bycoverage:
    input:
        expand(f"{RESULTDIR}/ngmlr/cutesv/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "DV>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """


rule lra_cutesv_bycoverage:
    input:
        expand(f"{RESULTDIR}/lra/cutesv/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "DV>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"CUTESV","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """


rule minimap2_sniffles_bycoverage:
    input:
        expand(f"{RESULTDIR}/minimap2/sniffles/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && cat <(cat ${{var}}| grep "^#") <(cat ${{var}}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{s}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """


rule ngmlr_sniffles_bycoverage:
    input:
        expand(f"{RESULTDIR}/ngmlr/sniffles/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && cat <(cat ${{var}}| grep "^#") <(cat ${{var}}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude}| sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{s}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """


rule lra_sniffles_bycoverage:
    input:
        expand(f"{RESULTDIR}/lra/sniffles/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && cat <(cat ${{var}}| grep "^#") <(cat ${{var}}| grep -vE "^#" | grep -v "0/0:" | grep -v -f {params.exclude} | sort -k1,1 -k2,2g) | bgzip -c > {params.vcfgz} && bcftools view -f PASS -i "DV>=${{s}}" -o {params.vcf} {params.vcfgz} && rm {params.vcfgz} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SNIFFLES","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule minimap2_svim_bycoverage:
    input:
        expand(f"{RESULTDIR}/minimap2/svim/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "SUPPORT>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | sed 's/DUP:TANDEM/DUP/g' | sed 's/DUP:INT/DUP/g'| bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule ngmlr_svim_bycoverage:
    input:
        expand(f"{RESULTDIR}/ngmlr/svim/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "SUPPORT>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | sed 's/DUP:TANDEM/DUP/g' | sed 's/DUP:INT/DUP/g'| bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule lra_svim_bycoverage:
    input:
        expand(f"{RESULTDIR}/lra/svim/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "SUPPORT>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | sed 's/DUP:TANDEM/DUP/g' | sed 's/DUP:INT/DUP/g'| bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"SVIM","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule pbmm2_pbsv_bycoverage:
    input:
        expand(f"{RESULTDIR}/pbmm2/pbsv/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && bcftools view -f PASS -i "AD[:1]>=${{s}}" ${{var}} | bcftools view -e 'GT[*]="RR"' | grep -v -f {params.exclude} | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"PBSV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """


rule minimap2_npinv_bycoverage:
    input:
        expand(f"{RESULTDIR}/minimap2/npinv/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/minimap2/npinv/npinv.minimap2.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/npinv_minimap2_prf1bycoverage.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/minimap2/npinv/tmp.vcf",
        vcfgz=f"{RESULTDIR}/minimap2/npinv/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/minimap2/npinv/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/minimap2/npinv",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && cat <(grep "^#" ${{var}}) <(grep -v -e "#" -e "_K" -e "HLA" -e "chrUn_" -e "_GL" -e "JH" -e "EB" -e "chrM" ${{var}}) | bcftools view -f PASS -i "DV>=${{s}}" | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"NPINV","minimap2"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """


rule ngmlr_npinv_bycoverage:
    input:
        expand(f"{RESULTDIR}/ngmlr/npinv/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/ngmlr/npinv/npinv.ngmlr.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/npinv_ngmlr_prf1bycoverage.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/ngmlr/npinv/tmp.vcf",
        vcfgz=f"{RESULTDIR}/ngmlr/npinv/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/ngmlr/npinv/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/ngmlr/npinv",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && cat <(grep "^#" ${{var}}) <(grep -v -e "#" -e "_K" -e "HLA" -e "chrUn_" -e "_GL" -e "JH" -e "EB" -e "chrM" ${{var}}) | bcftools view -f PASS -i "DV>=${{s}}" | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"NPINV","ngmlr"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """


rule lra_npinv_bycoverage:
    input:
        expand(f"{RESULTDIR}/lra/npinv/{{coverage}}/SI00001.vcf", coverage=["5X", "10X", "15X", "20X", "25X", "35X", "total"])
    output:
        f"{RESULTDIR}/lra/npinv/npinv.lra.bycoverage.tsv"  
    log:
        f"{LOGDIR}/results/npinv_lra_prf1bycoverage.log"
    threads: 1
    params:
        exclude = config["excludebed"],
        vcf=f"{RESULTDIR}/lra/npinv/tmp.vcf",
        vcfgz=f"{RESULTDIR}/lra/npinv/tmp.vcf.gz",
        vcfgztbi=f"{RESULTDIR}/lra/npinv/tmp.vcf.gz.tbi",
        outdir=f"{RESULTDIR}/lra/npinv",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tcoverage\ttool\taligner" > {output} && \
        cov="5:2 10:4 15:7 20:9 25:10 35:10 total:10" && files=$(ls {params.outdir}/*/SI00001.vcf | sort -V) && set ${{cov}} && \
        for var in ${{files}}; do c=$(echo ${{1}} | cut -f1 -d ":") && s=$(echo ${{1}} | cut -f2 -d ":") && cat <(grep "^#" ${{var}}) <(grep -v -e "#" -e "_K" -e "HLA" -e "chrUn_" -e "_GL" -e "JH" -e "EB" -e "chrM" ${{var}}) | bcftools view -f PASS -i "DV>=${{s}}" | bcftools view -e 'GT[*]="RR"' | bcftools sort > {params.vcf} && bgzip {params.vcf} && tabix {params.vcfgz} && truvari bench -f {params.ref} -b {params.truthvcf} -p 0 -c {params.vcfgz} -o {params.outdir}/truvari_coverage_${{c}} && cat {params.outdir}/truvari_coverage_${{c}}/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk -v su=${{c}} '{{OFS=FS="\t"}}{{print $0,su,"NPINV","lra"}}' >> {output} && rm {params.vcfgz} && rm {params.vcfgztbi} 2>>{log} && shift; done
        """

rule combine_results_bycoverage:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/{{tool}}/{{tool}}.{{aligner}}.bycoverage.tsv", tool=["cutesv", "svim", "sniffles", "npinv"], aligner=["minimap2", "ngmlr", "lra"]),
        f"{RESULTDIR}/pbmm2/pbsv/pbsv.pbmm2.bycoverage.tsv"
    output:
        f"{RESULTDIR}/SI00001.prf1.bycoverage.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/combine_prf1bycoverage.log"
    shell:
        """
        awk "FNR>1 || NR==1" {input} > {output} 2>{log}
        """

rule plot_results_bycoverage:
    input:
        f"{RESULTDIR}/SI00001.prf1.bycoverage.tsv"
    output:
        f"{RESULTDIR}/SI00001.prf1.bycoverage.pdf"
    threads: 1
    conda: "../envs/plot.yaml"    
    log:
        f"{LOGDIR}/results/plot_prf1bycoverage.log"
    shell:
        "Rscript {SCRIPTDIR}/prf1bycoverage.R {RESULTDIR} 2> {log}"


rule prf1_minimap2_cutesv_sniffles:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "sniffles"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_sniffles && cat {params.outdir}/truvari_minimap2_cutesv_sniffles/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_svim:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "svim"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_svim && cat {params.outdir}/truvari_minimap2_cutesv_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-svim","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_pbsv:    
    input:
        f"{RESULTDIR}/minimap2/cutesv/total/SI00001.clean.vcf",
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf",
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_pbsv && cat {params.outdir}/truvari_minimap2_cutesv_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_sniffles_pbsv:    
    input:
        f"{RESULTDIR}/minimap2/sniffles/total/SI00001.clean.vcf",
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf",
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_sniffles_pbsv && cat {params.outdir}/truvari_minimap2_sniffles_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"sniffles-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_svim_pbsv:
    input:
        f"{RESULTDIR}/minimap2/svim/total/SI00001.clean.vcf",
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf",
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_svim_pbsv && cat {params.outdir}/truvari_minimap2_svim_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"svim-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/SI00001.clean.vcf", tools=["sniffles", "svim"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_sniffles_svim && cat {params.outdir}/truvari_minimap2_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"sniffles-svim","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv","sniffles", "svim"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_sniffles_svim && cat {params.outdir}/truvari_minimap2_cutesv_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles-svim","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_sniffles_pbsv:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv","sniffles"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf"
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_sniffles_pbsv && cat {params.outdir}/truvari_minimap2_cutesv_sniffles_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_svim_pbsv:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv","svim"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf"
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_svim_pbsv && cat {params.outdir}/truvari_minimap2_cutesv_svim_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-svim-pbsv","minimap2"}}' >> {output}
        """ 


rule prf1_minimap2_sniffles_svim_pbsv:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/SI00001.clean.vcf", tools=["sniffles","svim"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf"
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_sniffles_svim_pbsv && cat {params.outdir}/truvari_minimap2_sniffles_svim_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"sniffles-svim-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_minimap2_cutesv_sniffles_svim_pbsv:    
    input:
        expand(f"{RESULTDIR}/minimap2/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv","sniffles","svim"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.clean.vcf"
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 4 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_minimap2_cutesv_sniffles_svim_pbsv && cat {params.outdir}/truvari_minimap2_cutesv_sniffles_svim_pbsv/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles-svim-pbsv","minimap2"}}' >> {output}
        """ 

rule prf1_ngmlr_cutesv_sniffles:    
    input:
        expand(f"{RESULTDIR}/ngmlr/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "sniffles"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_ngmlr_cutesv_sniffles && cat {params.outdir}/truvari_ngmlr_cutesv_sniffles/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles","ngmlr"}}' >> {output}
        """ 

rule prf1_ngmlr_cutesv_svim:    
    input:
        expand(f"{RESULTDIR}/ngmlr/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "svim"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_ngmlr_cutesv_svim && cat {params.outdir}/truvari_ngmlr_cutesv_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-svim","ngmlr"}}' >> {output}
        """ 

rule prf1_ngmlr_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/ngmlr/{{tools}}/total/SI00001.clean.vcf", tools=["sniffles", "svim"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_ngmlr_sniffles_svim && cat {params.outdir}/truvari_ngmlr_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"sniffles-svim","ngmlr"}}' >> {output}
        """ 

rule prf1_ngmlr_cutesv_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/ngmlr/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv","sniffles", "svim"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_ngmlr_cutesv_sniffles_svim && cat {params.outdir}/truvari_ngmlr_cutesv_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles-svim","ngmlr"}}' >> {output}
        """ 

rule prf1_lra_cutesv_sniffles:    
    input:
        expand(f"{RESULTDIR}/lra/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "sniffles"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_lra_cutesv_sniffles && cat {params.outdir}/truvari_lra_cutesv_sniffles/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles","lra"}}' >> {output}
        """ 

rule prf1_lra_cutesv_svim:    
    input:
        expand(f"{RESULTDIR}/lra/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv", "svim"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_lra_cutesv_svim && cat {params.outdir}/truvari_lra_cutesv_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-svim","lra"}}' >> {output}
        """ 

rule prf1_lra_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/lra/{{tools}}/total/SI00001.clean.vcf", tools=["sniffles", "svim"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 2 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_lra_sniffles_svim && cat {params.outdir}/truvari_lra_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"sniffles-svim","lra"}}' >> {output}
        """ 

rule prf1_lra_cutesv_sniffles_svim:    
    input:
        expand(f"{RESULTDIR}/lra/{{tools}}/total/SI00001.clean.vcf", tools=["cutesv","sniffles", "svim"]),
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
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        ls {input} > {params.tmp} && mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        SURVIVOR merge {params.tmp} 1000 3 1 1 0 50 {params.mergedvcf} && rm {params.tmp} && \
        bcftools sort -O z -o {params.mergedvcfgz} {params.mergedvcf} && bcftools index -t {params.mergedvcfgz} && rm {params.mergedvcf} && \
        truvari bench -f {params.ref} -b {params.truthvcf}  -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_lra_cutesv_sniffles_svim && cat {params.outdir}/truvari_lra_cutesv_sniffles_svim/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"cutesv-sniffles-svim","lra"}}' >> {output}
        """ 

rule prf1_minimap2_ngmlr_lra_consensus:    
    input:
        f"{RESULTDIR}/SI00001.consensus.vcf"
    output:
        f"{RESULTDIR}/calls_combined/consensus.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/consensus.log"
    params:
        outdir=f"{RESULTDIR}/calls_combined",
        tmp=f"{RESULTDIR}/calls_combined/minimap2_ngmlr_lra.txt",
        mergedvcfgz=f"{RESULTDIR}/calls_combined/minimap2_ngmlr_lra.vcf.gz",
        truthvcf=f"{VCFDIR}/hack.vcf.gz",
        ref=config["genome"]
    shell:
        """
        mkdir -p {params.outdir} && echo -e "precision\trecall\tf1\tprecision_gt\trecall_gt\tf1_gt\tsupport\ttool\taligner" > {output} && \
        bcftools sort -O z -o {params.mergedvcfgz} {input} && bcftools index -t {params.mergedvcfgz} && \
        truvari bench -f {params.ref} -b {params.truthvcf} --passonly -p 0 -c {params.mergedvcfgz} -o {params.outdir}/truvari_consensus && cat {params.outdir}/truvari_consensus/summary.txt | tail -n+2 | head -n -1 | sed "s/^[ \t]*//" | grep -E "precision|recall|f1|gt_precision|gt_recall|gt_f1" | cut -d ":" -f 2 | sed "s/,//g" | sed "s/^ //g" | paste -s -d "\t" | awk '{{OFS=FS="\t"}}{{print $0,10,"consensus","consensus"}}' >> {output}
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
        f"{RESULTDIR}/calls_combined/minimap2_cutesv_sniffles_svim_pbsv.tsv",
        f"{RESULTDIR}/calls_combined/consensus.tsv"
    output:
        f"{RESULTDIR}/SI00001.prf1.bycombo.tsv"
    threads: 1
    log:
        f"{LOGDIR}/results/combine_prf1bycombo.log"
    shell:
        """
        awk "FNR>1 || NR==1" {input} > {output} 2>{log}
        """

rule plot_results_bycombo:
    input:
        f"{RESULTDIR}/SI00001.prf1.bycombo.tsv"
    output:
        f"{RESULTDIR}/SI00001.prf1.bycombo.pdf"
    threads: 1
    conda: "../envs/plot.yaml"    
    log:
        f"{LOGDIR}/results/plot_prf1bycombo.log"
    shell:
        "Rscript {SCRIPTDIR}/prf1bycombo.R {RESULTDIR} 2> {log}"

rule combine_stats_by_length:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/{{tool}}/total/{{tool}}.{{aligner}}.bysupport.tsv", tool=["cutesv", "svim", "sniffles", "npinv"], aligner=["minimap2", "ngmlr","lra"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/pbsv.pbmm2.bysupport.tsv"
    output:
        f"{RESULTDIR}/SI00001.tpfpfn.bylength.tsv"
    log:
        f"{LOGDIR}/results/combine_stats_by_length.log"
    shell:
        """
        echo -e "size\tsvtype\tvartype\ttool\taligner" > {output} && \
        find {RESULTDIR} -wholename '*truvari_supporting_10/svlen.info.tsv*' -exec cat {{}} + >> {output} 2> {log}
        """

rule bnds_prf1:
    input:
        expand(f"{RESULTDIR}/{{aligner}}/{{tools}}/total/SI00001.bnd.bedpe", aligner=["minimap2", "ngmlr", "lra"], tools=["sniffles", "svim", "cutesv"]),
        f"{RESULTDIR}/pbmm2/pbsv/total/SI00001.bnd.bedpe"
    output:
        f"{RESULTDIR}/SI00001.bnds.prf1.tsv"
    threads: 1
    conda: "../envs/pybedtools.yaml"
    log:
        f"{LOGDIR}/results/bnds_prf1"
    params:
        truthbedpe=f"{VCFDIR}/SI00001.bnd.bedpe"
    shell:
        """
        for var in {input}; do t=$(basename $(dirname $(dirname $(readlink -f ${{var}})))) && a=$(basename $(dirname $(dirname $(dirname $(readlink -f ${{var}}))))) && python {SCRIPTDIR}/evalbnd.py {params.truthbedpe} ${{var}} | awk -v al=${{a}} -v to=${{t}} '{{FS=OFS="\t"}}''{{print $0,to,al}}' >> {output} 2>>{log}; done && sed -i "s/pbmm2/minimap2/g" {output}
        """


rule plot_stats_by_length:
    input:
        f"{RESULTDIR}/SI00001.tpfpfn.bylength.tsv",
        f"{RESULTDIR}/SI00001.bnds.prf1.tsv"
    output:
        f"{RESULTDIR}/SI00001.tpfpfn.bylength.pdf",
        f"{RESULTDIR}/SI00001.svtype.bylength.pdf",
        f"{RESULTDIR}/SI00001.prf1.bysvtype.pdf"
    threads: 1
    conda: "../envs/plot.yaml"    
    log:
        f"{LOGDIR}/results/plot_stats_by_length"
    shell:
        "Rscript {SCRIPTDIR}/tpfpfnbysize.R {RESULTDIR} 2> {log}"

