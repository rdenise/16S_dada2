# Snakemake workflow: Virome pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.24.2-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Aim

This workflow uses DADA2 to process RNA 16S reads. The pipeline consists of three main rules:

quality: This rule takes the raw reads and performs quality control, trimming, and filtering.
all: This rule takes the preprocessed reads and performs denoising, chimera removal, and generates the final feature table and representative sequences.

## Dependencies

Snakemake and Snakedeploy are best installed via the [Mamba package manager](https://github.com/mamba-org/mamba) (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). For other options see [here](https://github.com/mamba-org/mamba).

To install Mamba by conda, run

```shell
conda install mamba -n base -c conda-forge
```

Given that Mamba is installed, run 

```shell
mamba create -c bioconda -c conda-forge --name snakemake snakemake snakedeploy
```

to install both Snakemake and Snakedeploy in an isolated environment. 


## Configuring the pipeline

Before running the pipeline, you need to modify the config.yaml file to match your specific project requirements.

reads_folder: This is the path to the folder containing your raw reads.
reads_identifier: This is the string that identifies the pair end name in the raw read files.
output_folder: This is the path to the folder where you want to store the results.
silva_nr_train_set: This is the path to the SILVA NR99 reference database.
lenght_reads: Th wanted size to trim the reads

## Running the pipeline

To run the pipeline, navigate to the directory containing the Snakefile and type the following command:

### First look at the quality of your reads using 

```bash
snakemake quality --use-conda --cores 10
```

Using this file you will identify the size you reads sould be trimmed (normally it is around 200-250bp). The quality of the sequence should be above 20.

### Once you identify that the size to trim

Change in the config file the size (important remember that the untrimmed reads contains the adaptor at the begining so think about the size without the adaptor)

```bash
snakemake --use-conda --cores 10
```

This command will automatically download and install the necessary dependencies using Conda, and execute the workflow.

## Output

The pipeline will generate several output files, including:

Quality control reports
Trimmed and filtered read files
DADA2 output files, including feature table and representative sequences

The final output files will be located in the `output_folder` specified in the config.yaml file.

