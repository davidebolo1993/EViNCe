# NA24385

Oxford Nanopore PromethION ultra long reads of the Ashkenazim son [HG002/NA24385](ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion). GIAB variant calls for the same individual are the [ground truth](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/)


## Get Data

``` bash
bash scripts/getdata.sh
```

## Alignment

``` bash
snakemake --cores 8 alignments/GM24385.{minimap2,ngmlr}.srt.bam.bai
```