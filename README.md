# GLIMPSE Example Snakemake Workflow

## Overview

Most of this if from old glimpse pipelien that exists here:https://github.com/odelaneau/GLIMPSE. The current script was modified for hg38 and to include respective reference panel. 

This workflow provides an end-to-end example of an imputation with GLIMPSE. The only dependency required is snakemake itself (`conda create -c conda-forge -c bioconda -n snakemake snakemake`). All other tools (bwa, gatk, rasusa, samtools, bcftools, and GLIMPSE) are integrated into the workflow using conda and singularity. The workflow performs the following:

1. Downloads the human reference genome (GRCh37).
2. Indexes the genome to create files needed downstream (samtools and bwa).
3. Downloads example reads from NA12878 (high coverage sequencing data from https://www.ebi.ac.uk/ena/browser/view/PRJEB8596).
4. Subsamples NA12878 to 30x, 1x, and 0.1x coverage.
5. For each coverage, map (bwa), sort, mark dupliates.
6. Call variants with GATK for the 30x sample only (useful for optional comparison purposes downstream, not part of this pipeline).
7. Call variants with BCFtools as described in the GLIMPSE documentation for downstream GLIMPSE imputation.
8. Create the GLIMPSE imputation panel from 1000 Genomes data.
9. Perform low-coverage WGS with GLIMPSE.

![](dag.svg)

## Usage

This workflow must be run in two steps, as GLIMPSE requires "chunk" files which must be precomputed prior to executing the main workflow. 

Before starting the main workflow, you will need to generate the chunkfiles from the 1000 genomes sites.vcf.gz files. First execute this rule, which will download the necessary reference data and create chunks that the imputation workflow uses later. These must be pre-populated and cannot be made concurrently with the main workflow.

```sh
snakemake -j10 --use-singularity --use-conda -s workflow/Make_chunks.smk
```

Run the entire workflow after this step completes:

```sh
snakemake -j20 --use-singularity --use-conda
```

