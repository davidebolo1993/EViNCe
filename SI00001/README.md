# SI00001


### Create Data

We first need to install [VISOR](https://github.com/davidebolo1993/VISOR). Instructions are provided [here](https://davidebolo1993.github.io/visordoc/installation/installation.html#installation-from-source).

``` bash
cd workflow/scripts
bash getdata.sh
cd -
```

### Align FASTQ to GRCh38

``` bash
snakemake align --cores 20 
```

### Calculate/plot coverage statistics

``` bash
snakemake depth --cores 20
```

### Subsample original alignments

``` bash
snakemake subsample --cores 20
```

