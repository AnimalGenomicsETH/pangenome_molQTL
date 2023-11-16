[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![snakemaker](https://github.com/AnimalGenomicsETH/pangenome_molQTL/actions/workflows/snakemaker.yaml/badge.svg?branch=main)](https://github.com/AnimalGenomicsETH/pangenome_molQTL/actions/workflows/snakemaker.yaml)
# Pangenome genotyping of structural variants and molecular QTL mapping

Structural variants are known to play a large role in expression and splicing QTL.
However, confidently calling stuctural variants in a sufficiently large population for association mapping is hard.
Here we use [PanGenie](https://github.com/eblerjana/pangenie) to genotype a larger cohort (100s) of short reads using an accurate pangenome panel from haplotype-resolved assemblies (10s).

### Usage

There are broadly three phases
  - Pangenome panel (creating & genotyping)
  - Variant analysis (statistics, linkage disequibrium, SV overlap, etc.) 
  - Association mapping of e/sQTL


An example of the input needed is given in the `config/example.yaml`, broadly requiring
  - haplotype-resolved assemblies for pangenome panel creation
  - small variants to supplement pangenome panel
  - any HiFi samples to test SV completeness
  - gene expression/splicing files and covariates for molecular QTL mapping

Running with
```
snakemake --configfile config/example.yaml
```
Will execute the following DAG

![workflow](https://github.com/AnimalGenomicsETH/pangenome_molQTL/assets/29678761/bb0c73ca-fc31-4319-95e2-485da93f655a)

producing the major output files (e.g., accuracy comparison of PanGenie vs DeepVariant, SV overlap with Jasmine, conditional QTL analysis with QTLtools, etc.), which can then be independently analysed further.
Many of these steps are computationally intensive, especially with many samples to genotype, and so effectively require some form of HPC.

### Citation

The preprint associated with this work can be found [here](https://www.biorxiv.org/content/10.1101/2023.06.21.545879v1).
> **Pangenome genotyped structural variation improves molecular phenotype mapping in cattle**
> 
> Alexander S. Leonard, Xena M. Mapel, Hubert Pausch

### Note
Many of the parameters are tuned to run for our data and on the ETH Euler cluster, using for example a forked version of the LSF snakemake [profile](https://github.com/AnimalGenomicsETH/snakemake_lsf), so it may take some modifying to work smoothly in different contexts. Many tools are assumed to be available in $PATH, but all are freely available.

