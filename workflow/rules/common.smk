##########################################################################
##########################################################################
##
##                                Library
##
##########################################################################
##########################################################################

import os, sys
import pandas as pd
import numpy as np
from snakemake.utils import validate
import glob

##########################################################################
##########################################################################
##
##                               Functions
##
##########################################################################
##########################################################################


def get_final_output(outdir):
    """
    Generate final output name
    """
    final_output = []

    final_output += (
        os.path.join(
            outdir,
            "reports",
            "qc",
            "multiqc_report.trimmed.html",
        ),
    )

    final_output += (
        os.path.join(
            outdir,
            "results",
            "count_table_from_dada2.csv"
        ),
    )

    return final_output


##########################################################################


def create_folder(mypath):
    """
    Created the folder that I need to store my result if it doesn"t exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return


##########################################################################
##########################################################################
##
##                                Variables
##
##########################################################################
##########################################################################

# Validation of the config.yaml file
validate(config, schema="../schemas/config.schema.yaml")

##########################################################################
##########################################################################
##
##                        Core configuration
##
##########################################################################
##########################################################################

## Store some workflow metadata
config["__workflow_basedir__"] = workflow.basedir
config["__workflow_basedir_short__"] = os.path.basename(workflow.basedir)
config["__workflow_workdir__"] = os.getcwd()

if workflow.config_args:
    tmp_config_arg = '" '.join(workflow.config_args).replace("=", '="')
    config["__config_args__"] = f' -C {tmp_config_arg}"'
else:
    config["__config_args__"] = ""

with open(os.path.join(workflow.basedir, "../config/VERSION"), "rt") as version:
    url = "https://github.com/rdenise/viral_detection/releases/tag"
    config["__workflow_version__"] = version.readline()
    config["__workflow_version_link__"] = f"{url}/{config['__workflow_version__']}"


##########################################################################
##########################################################################
##
##                           Options
##
##########################################################################
##########################################################################

# Result folder
OUTPUT_FOLDER = config["output_folder"]
# Adding to config for report
config["__output_folder__"] = os.path.abspath(OUTPUT_FOLDER)

# Get fastq folder
FASTQ_FOLDER = config["reads_folder"]

# Get fastq separator
FASTQ_SEP = config["reads_identifier"]

# Get the samples
(FASTQ_SAMPLES, FASTQ_EXT,) = glob_wildcards(
    os.path.join(
        FASTQ_FOLDER, 
        "{fastq_files}" + FASTQ_SEP + "1" + "{fastq_extension}"
    )
)

# Only get one extension
FASTQ_EXT = FASTQ_EXT[0] if "md5" not in FASTQ_EXT[0] else FASTQ_EXT[0][:-4]

# print(get_final_output(
#             outdir=OUTPUT_FOLDER,
#             reads=reads
#         ))
