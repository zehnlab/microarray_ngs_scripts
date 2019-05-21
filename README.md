This repository provides the scripts of Microarray & NGS data analysis for the paper: "Alfei et al. Tox reinforces the phenotype and longevity of dysfunctional T cell populations during chronic viral infection."

**code**: microarray analysis + gene differentiation expression analysis + heatmap visualization;

**data**: readcount tables + design tables;

**output**: differential expressed gene lists + heatmaps.

To create the outputs from the data and the code you can use snakemake as follows:

## Usage

### Step 1: Snakemake installation

Please follow the [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

### Step 2: Clone the repo

`git clone https://github.com/zehnlab/microarray_ngs_scripts`

### Step 3: Run the workflow

go the cloned dir:
`cd microarray_ngs_scripts`

To generate the RNA differential expression run: 
`snakemake --use-conda DGE_RNAseq`

To generate the heatmaps run: 
`snakemake --use-conda generate_heatmaps`