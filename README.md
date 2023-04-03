# Snakemake workflow: 16S RNA workflow with dada2

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.24.2-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/619899863.svg)](https://zenodo.org/badge/latestdoi/619899863)


## Aim

This workflow uses DADA2 to process RNA 16S reads. The pipeline consists of three main rules:

quality: This rule takes the raw reads and performs quality control, trimming, and filtering.  
all: This rule takes the preprocessed reads and performs denoising, chimera removal, and generates the final feature table and representative sequences.  

This workflow is based on the DADA2 Rscript made by Thomaz Bastiaanssen (UCC, APC Microbiome). A copy of the script is included in the scripts folder. A guidebook on how to calculate alpha and beta diversity is available at this url https://github.com/thomazbastiaanssen/Tjazi.

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

### Deploy workflow

Given that Snakemake and Snakedeploy are installed and available (see Step 1), the workflow can be deployed as follows.

First, create an appropriate project working directory on your system in the place of your choice as follow (note to change the path and file name to the one you want to create): : 

```shell
mkdir path/to/project-workdir
```

Then go to your project working directory as follow:

```shell
cd path/to/project-workdir
```

In all following steps, we will assume that you are inside of that directory.

Second, run 

```shell
snakedeploy deploy-workflow https://github.com/rdenise/detection_virus_metagenomes . --tag 0.0.1
```

Snakedeploy will create two folders `workflow` and `config`. The former contains the deployment of the chosen workflow as a [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows), the latter contains configuration files which will be modified in the next step in order to configure the workflow to your needs. Later, when executing the workflow, Snakemake will automatically find the main `Snakefile` in the `workflow` subfolder.


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

A folder containing your work will be created:

```
[output_folder]                            <- Main results folder
├── databases                              <- Folder containing databases used in the analysis
│   └── reads_trimmed                      <- Folder containing trimmed read FASTQ files
├── logs                                   <- Folder containing log files for each analysis step
├── processed_files                        <- Folder containing processed files resulting from the analysis
│   └── dada2                              <- Folder containing output files from DADA2 analysis pipeline
│       ├── denoised                       <- Folder containing denoised sequences
│       ├── filtered-pe                    <- Folder containing filtered paired-end reads
│       ├── learn-errors                   <- Folder containing learned error rates
│       ├── merged                         <- Folder containing merged paired-end reads
│       └── uniques                        <- Folder containing unique sequences
├── reports                                <- Folder containing analysis reports
│   ├── dada2                              <- Folder containing DADA2 analysis pipeline reports
│   │   ├── filter-trim-pe                 <- Report on quality filtering and trimming of paired-end reads
│   │   └── quality-profile                <- Report on sequence quality profile
│   └── qc                                 <- Folder containing quality control reports
│       ├── fastqc                         <- FastQC report files for each input read file
│       ├── multiqc_report.trimmed_data    <- MultiQC report for trimmed reads
│       └── multiqc_report.untrimmed_data  <- MultiQC report for untrimmed reads
└── results                                <- Folder containing final analysis results

```