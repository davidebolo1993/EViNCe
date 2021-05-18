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

### Call SVs with cuteSV

``` bash
snakemake cutesv --cores 20
```

### Call SVs with sniffles

``` bash
snakemake sniffles --cores 20
```

### Call SVs with svim

``` bash
snakemake svim --cores 20
```

### Call INVs with npinv

``` bash
snakemake npinv --cores 20
```

### Compute/plot upset plots and plot variants per chromosome

``` bash
snakemake sumsvs --cores 20
```

