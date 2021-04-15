# EViNCe
EValuation of Nanopore (SV) Callers

The aim is to evaluate the performances of several long-read structural variant (SV) callers on long reads from Oxford Nanopore Technologies. This will include tests both on real and simulated datasets.

## SV callers

To date, the following methods (versions) are tested:

- [svim](https://github.com/eldariont/svim) (1.4.2)
- [sniffles](https://github.com/fritzsedlazeck/Sniffles) (1.0.12)
- [pbsv](https://github.com/PacificBiosciences/pbsv) (2.4.0)
- [cuteSV](https://github.com/tjiangHIT/cuteSV) (1.0.10)
- [npInv](https://github.com/haojingshao/npInv) (1.24)

Long-read aligners tested:

- [minimap2](https://github.com/lh3/minimap2) (2.17)
- [ngmlr](https://github.com/philres/ngmlr) (0.2.7)
- [lra](https://github.com/ChaissonLab/LRA) (1.1.2)
- [pbmm2](https://github.com/PacificBiosciences/pbmm2) (1.3.0) 

Other tools used by the pipeline:

- [samtools](https://github.com/samtools/samtools) (1.9), using [htslib](https://github.com/samtools/htslib) (1.9)
- [bcftools](https://github.com/samtools/bcftools) (1.9), using [htslib](https://github.com/samtools/htslib) (1.9)
- [truvari](https://github.com/spiralgenetics/truvari) (2.0.0-dev)
- [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR.git) (1.0.7)
- [VISOR](https://github.com/davidebolo1993/VISOR) (1.1)

Pipelines are built using [Snakemake](https://snakemake.readthedocs.io/en/stable) to improve readability/reproducibility of the results. This also allows to install most of the tools required for each command on-the-fly, using the --use-conda flag. 
