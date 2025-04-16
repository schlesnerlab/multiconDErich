# multiconDErich

[![Tests](https://github.com/schlesnerlab/multicondition-deseq2-enrichment/actions/workflows/main.yaml/badge.svg)](https://github.com/schlesnerlab/multicondition-deseq2-enrichment/actions/workflows/main.yaml)

Snakemake workflow for running differential expression and enrichment analyses for experimental 
setups with multiple groups or multiple conditions. Results are provided as a set of HTML files
with plots and result description both for each comparison defined in the config file. 

Workflow was derived from https://github.com/snakemake-workflows/rna-seq-star-deseq2, with 
alignment steps removed. 


## Installation

### Required dependencies

The quickest way to get up and running with this workflow is using
- [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
	- [mamba](https://github.com/mamba-org/mamba)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)
- [pandas]()

### Installation


```bash
git clone https://github.com/schlesnerlab/multicondition-deseq2-enrichment
cd multicondition-deseq-enrichment
```

You can install snakemake with the yaml provided in `workflow/snakemake.yaml`
with 

```bash
conda env create -n snakemake -f workflow/snakemake.yaml
conda activate snakemake
```

In case you use snakedeploy you could also deploy this workflow with
```bash
snakedeploy deploy-workflow https://github.com/schlesnerlab/multicondition-deseq2-enrichment multicondition-deseq2-enrichment --branch main 
``` 

### Carnival usage 

If you wish to use analysis provided by Carnival you need to provide a solver
supported by CARNIVAL. Right now the following solvers are supported:

- cplex

For academic users the cplex can be downloaded [here](https://www.ibm.com/products/ilog-cplex-optimization-studio). Then the path to the cplex executable needs to 
added to the config file. 

## Usage

### Setting up config and sample sheet

To use this workflow you will need:
- RNAseq quantified data (f.e. counts) as a tsv or csv file
- a config yaml file for your project to control the workflow
- a tsv file with the sample information for the samples in the count file.

Further details on the configuration can be found in the config README


### Supported input data

- gene count matrix
