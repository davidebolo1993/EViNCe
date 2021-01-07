# EViNCe
EValuation of Nanopore (SV) Callers

This is a work in progress. The aim is to evaluate the performances of several long-read structural variant (SV) callers on long reads from Oxford Nanopore Technologies. This will include tests both on real and simulated datasets.

## Tools

For the time being, the following methods (versions) will be tested. Maybe others will be added in the future.

- [svim](https://github.com/eldariont/svim) (1.4.2)
- [sniffles](https://github.com/fritzsedlazeck/Sniffles) (1.0.12)
- [nanosv](https://github.com/mroosmalen/nanosv) (1.2.4)
- [cuteSV](https://github.com/tjiangHIT/cuteSV) (1.0.10)
- [nanomonsv](https://github.com/friend1ws/nanomonsv) (0.2.0)

Tools used in the analysis:

- [minimap2](https://github.com/lh3/minimap2) (2.17)
- [ngmlr](https://github.com/philres/ngmlr) (0.2.7)
- ~~[samtools](https://github.com/samtools/samtools) (1.11), using [htslib](https://github.com/samtools/htslib) (1.11)~~ , due to "incomplete aux fields" errors occuring when sorting alignments from ngmlr with samtools 1.10/1.11 (there are several open issues about this, https://github.com/philres/ngmlr/issues/89 and https://github.com/philres/ngmlr/issues/86 for instance). Downgrading seems to work.
- [samtools](https://github.com/samtools/samtools) (1.9), using [htslib](https://github.com/samtools/htslib) (1.9)
- ~~[bcftools](https://github.com/samtools/bcftools) (1.11), using [htslib](https://github.com/samtools/htslib) (1.11)~~
- [bcftools](https://github.com/samtools/bcftools) (1.9), using [htslib](https://github.com/samtools/htslib) (1.9)
- [truvari](https://github.com/spiralgenetics/truvari) (2.0.0-dev)
- [VISOR](https://github.com/davidebolo1993/VISOR) (1.1)

Pipelines are built using [Snakemake](https://snakemake.readthedocs.io/en/stable) to improve readability/reproducibility of the results. Suggestions are very much 
appreciated.
