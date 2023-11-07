[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![snakemaker](https://github.com/AnimalGenomicsETH/pangenome_molQTL/actions/workflows/snakemaker.yaml/badge.svg?branch=main)](https://github.com/AnimalGenomicsETH/pangenome_molQTL/actions/workflows/snakemaker.yaml)
# Pangenome genotyping of structural variants and molecular QTL mapping

Structural variants are known to play a large role in expression and splicing QTL.
However, confidently calling stuctural variants in a sufficiently large population for association mapping is hard.
Here we use [PanGenie](https://github.com/eblerjana/pangenie) to genotype a larger cohort (100s) of short reads using an accurate pangenome panel from haplotype-resolved assemblies (10s).

### Usage

Three phases
  - Pangenome panel (creating & genotyping)
  - Variant analysis (statistics, linkage disequibrium, SV overlap, etc.) 
  - Association mapping of e/sQTL

Runs with
```
snakemake --configfile config/config.yaml
```

### Citation

The preprint associated with this work can be found [here](https://www.biorxiv.org/content/10.1101/2023.06.21.545879v1).
> **Pangenome genotyped structural variation improves molecular phenotype mapping in cattle**
> 
> Alexander S. Leonard, Xena M. Mapel, Hubert Pausch

### Note
Many of the parameters are tuned to run for our data and on the ETH Euler cluster, using for example a forked version of the LSF snakemake [profile](https://github.com/AnimalGenomicsETH/snakemake_lsf), so it may take some modifying to work smoothly in different contexts. Many tools are assumed to be available in $PATH, but all are freely available.
