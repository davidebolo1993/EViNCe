# EViNCe
EValuation of Nanopore (SV) Callers

This is a work in progress. The aim is to evaluate the performances of several long-read structural variant (SV) callers on long reads from Oxford Nanopore Technologies. This will include tests both on real and simulated datasets.

## Tools

For the time being, the following methods (versions) will be tested. Maybe others will be added in the future.

- [svim](https://github.com/eldariont/svim) (1.4.2)
- [sniffles](https://github.com/fritzsedlazeck/Sniffles) (1.0.12)
- [nanosv](https://github.com/mroosmalen/nanosv) (1.2.4)
- [smoove](https://github.com/brentp/smoove) (0.2.6), which includes best practices for [lumpy](https://github.com/arq5x/lumpy-sv).
- [cuteSV](https://github.com/tjiangHIT/cuteSV) (1.0.9)
- [nanovar](https://github.com/cytham/nanovar) (1.3.8)
- [delly](https://github.com/dellytools/delly) (0.8.6), which has recently introduced a method for long-read calling

In addition to those above, other tools are required in PATH:

- [minimap2](https://github.com/lh3/minimap2) (1.11)
- [ngmlr](https://github.com/philres/ngmlr) (0.2.7)
- ~~[samtools](https://github.com/samtools/samtools) (1.11), using [htslib](https://github.com/samtools/htslib) (1.11)~~ ,due to "incomplete aux fields" errors occuring when sorting alignments from ngmlr with samtools 1.10/1.11 (there are several open issues about this, https://github.com/philres/ngmlr/issues/89 and https://github.com/philres/ngmlr/issues/86 for instance). Downgrading seems to work.
- [samtools](https://github.com/samtools/samtools) (1.9), using [htslib](https://github.com/samtools/htslib) (1.9)
- ~~[[bcftools](https://github.com/samtools/bcftools) (1.11), using [htslib](https://github.com/samtools/htslib) (1.11)~~
- [[bcftools](https://github.com/samtools/bcftools) (1.9), using [htslib](https://github.com/samtools/htslib) (1.9)
- [truvari](https://github.com/spiralgenetics/truvari) (2.0.0-dev)
- [VISOR](https://github.com/davidebolo1993/VISOR) (1.1)

I'm gonna try to use [Snakemake](https://snakemake.readthedocs.io/en/stable) as much as possible to improve readability/reproducibility of the pipeline I use, even if I'm more confortable at scripting in bash. Suggestions are very much appreciated.
