# NA24385

Oxford Nanopore PromethION ultra long reads of the Ashkenazim son [HG002/NA24385](https://ftp.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/). GIAB variant calls for the same individual are also [available](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/) and can be used as ground truth.


## Download Data

``` bash
cd workflow/scripts
bash getdata.sh
cd -
```

## Align FASTQ to hs37d5

``` bash
snakemake align --cores 20 
```

## Calculate/plot coverage statistics

``` bash
snakemake depth --cores 20
```

## Subsample original alignments

``` bash
snakemake subsample --cores 20
```

## Call SVs with cuteSV

``` bash
snakemake cutesv --cores 20
```

## Call SVs with sniffles

``` bash
snakemake sniffles --cores 20
```

## Call SVs with svim

``` bash
snakemake svim --cores 20
```

