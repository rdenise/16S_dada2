# Make sure that you set the `truncLen=` option in the rule `dada2_filter_and_trim_pe` according
# to the results of the quality profile checks (after rule `dada2_quality_profile_pe` has finished on all samples).
# If in doubt, check https://benjjneb.github.io/dada2/tutorial.html#inspect-read-quality-profiles
##########################################################################
##########################################################################

rule dada2_quality_profile_pe:
    input:
        # FASTQ file without primer sequences
        expand(os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "reads_trimmed",
                "paired",
                "{{sample}}_R{orientation}.fastq.gz",
            ),
            orientation=[1,2]
        )
    output:
        os.path.join(
                OUTPUT_FOLDER,
                "reports",
                "dada2",
                "quality-profile",
                "{sample}-quality-profile.png"
        )
    log:
        os.path.join(
                OUTPUT_FOLDER,
                "logs",
                "dada2",
                "quality-profile",
                "{sample}-quality-profile.log"
        )
    wrapper:
        "v1.22.0/bio/dada2/quality-profile/wrapper.R"

##########################################################################
##########################################################################


rule dada2_filter_trim_pe:
    input:
        # Paired-end files without primer sequences
        fwd=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "reads_trimmed",
            "paired",
            "{sample}_R1.fastq.gz",
        ),
        rev=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "reads_trimmed",
            "paired",
            "{sample}_R2.fastq.gz",
        )
    output:
        filt=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "dada2",
            "filtered-pe",
            "{sample}_R1.fastq.gz"
        ),
        filt_rev=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "dada2",
            "filtered-pe",
            "{sample}_R2.fastq.gz"
        ),
        stats=os.path.join(
                OUTPUT_FOLDER,
                "reports",
                "dada2",
                "filter-trim-pe",
                "{sample}.tsv"
        )
    params:
        # Set the maximum expected errors tolerated in filtered reads
        maxEE=2,
        # Set the number of kept bases in forward and reverse reads
        truncLen=config["dada2"]["truncLen"],
        truncQ=2,
        maxN=0,
        verbose=True,
    log:
        os.path.join(
                OUTPUT_FOLDER,
                "logs",
                "dada2",
                "filter-trim-pe",
                "{sample}-filter-trim-pe.log"
        )
    threads: 5 # set desired number of threads here
    wrapper:
        "v1.22.0/bio/dada2/filter-trim/wrapper.R"

##########################################################################
##########################################################################


rule dada2_learn_errors:
    input:
        # Quality filtered and trimmed forward FASTQ files (potentially compressed)
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "filtered-pe",
                "{sample}_{{orientation}}.fastq.gz"
            ),
            sample=FASTQ_SAMPLES,
        )
    output:
        err=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "learn-errors",
                "model_{orientation}.RDS"
            ),# save the error model
        plot=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "learn-errors",
                "errors_{orientation}.png"
            ),# plot observed and estimated rates
    params:
        nbases=1e8
    log:
        os.path.join(
                OUTPUT_FOLDER,
                "logs",
                "dada2",
                "learn-errors",
                "learn-errors_{orientation}.log"
        )
    threads: 5 # set desired number of threads here
    wrapper:
        "v1.19.0/bio/dada2/learn-errors/wrapper.R"

##########################################################################
##########################################################################

rule dada2_dereplicate_fastq:
    input:
        # Quality filtered FASTQ file
        os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "filtered-pe",
                "{fastq}.fastq.gz"
        )
    output:
        # Dereplicated sequences stored as `derep-class` object in a RDS file
        os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "uniques",
                "{fastq}.RDS"
        )
    log:
        os.path.join(
                OUTPUT_FOLDER,
                "logs",
                "dada2",
                "dereplicate-fastq",
                "{fastq}.log"
        )
    wrapper:
        "v1.22.0/bio/dada2/dereplicate-fastq/wrapper.R"


##########################################################################
##########################################################################

rule dada2_sample_inference:
    input:
        # Dereplicated (aka unique) sequences of the sample
        derep=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "uniques",
                "{sample}_{orientation}.RDS"
        ),
        err=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "learn-errors",
                "model_{orientation}.RDS"
        ) # Error model
    output:
        os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "denoised",
                "{sample}_{orientation}.RDS"
        ) # Inferred sample composition
    log:
        os.path.join(
                OUTPUT_FOLDER,
                "logs",
                "dada2",
                "sample-inference",
                "{sample}_{orientation}.log"
        )
    threads: 5 # set desired number of threads here
    wrapper:
        "v1.22.0/bio/dada2/sample-inference/wrapper.R"


##########################################################################
##########################################################################

rule dada2_merge_pairs:
    input:
        dadaF=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "denoised",
                "{sample}_R1.RDS"
        ),# Inferred composition
        dadaR=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "denoised",
                "{sample}_R2.RDS"
        ),
        derepF=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "uniques",
                "{sample}_R1.RDS"
        ),# Dereplicated sequences
        derepR=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "uniques",
                "{sample}_R2.RDS"
        )
    output:
        os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "merged",
                "{sample}.RDS"
        ) # Merged composition
    log:
        os.path.join(
                OUTPUT_FOLDER,
                "logs",
                "dada2",
                "merge_pairs",
                "{sample}.log"
        )
    threads: 5 # set desired number of threads here
    wrapper:
        "v1.19.0/bio/dada2/merge-pairs/wrapper.R"

##########################################################################
##########################################################################

rule dada2_make_table_pe:
    input:
        # Merged composition
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "merged",
                "{sample}.RDS"
            ), 
            sample=FASTQ_SAMPLES,
        )
    output:
        os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "seqTab-pe.RDS"
        ) # Sequence table
    params:
        names=FASTQ_SAMPLES, # Sample names instead of paths
        orderBy="nsamples" # Change the ordering of samples
    log:
        os.path.join(
                OUTPUT_FOLDER,
                "logs",
                "dada2",
                "make-table",
                "make-table-pe.log"
        )
    threads: 5 # set desired number of threads here
    wrapper:
        "v1.19.0/bio/dada2/make-table/wrapper.R"

##########################################################################
##########################################################################

rule dada2_remove_chimeras:
    input:
        os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "seqTab-pe.RDS"
        ) # Sequence table
    output:
        os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "seqTab.nochimeras.RDS"
        ) # Chimera-free sequence table
    params:
        method="consensus",
    log:
        os.path.join(
                OUTPUT_FOLDER,
                "logs",
                "dada2",
                "remove-chimeras",
                "remove-chimeras.log"
        )
    threads: 5 # set desired number of threads here
    wrapper:
        "v1.19.0/bio/dada2/remove-chimeras/wrapper.R"

##########################################################################
##########################################################################

rule dada2_collapse_nomismatch:
    input:
        os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "seqTab.nochimeras.RDS"
        ) # Chimera-free sequence table
    output:
        os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "seqTab.collapsed.RDS"
        ) # Collapsed sequence table
    log:
        os.path.join(
                OUTPUT_FOLDER,
                "logs",
                "dada2",
                "collapse_nomismatch",
                "collapse_nomismatch.log"
        )
    threads: 5 # set desired number of threads here
    wrapper:
        "v1.19.0/bio/dada2/collapse-nomismatch/wrapper.R"

##########################################################################
##########################################################################

rule dada2_assign_taxonomy:
    input:
        seqs=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "seqTab.nochimeras.RDS"
        ), # Chimera-free sequence table
        refFasta=config["silva_nr_train_set"] # Reference FASTA for taxonomy
    output:
        os.path.join(
                OUTPUT_FOLDER,
                "results",
                "taxa.RDS"
        ) # Taxonomic assignments
    log:
        os.path.join(
                OUTPUT_FOLDER,
                "logs",
                "dada2",
                "assign_taxonomy",
                "assign_taxonomy.log"
        )
    threads: 5 # set desired number of threads here
    wrapper:
        "v1.19.0/bio/dada2/assign-taxonomy/wrapper.R"

##########################################################################
##########################################################################

rule results_table_format:
    input:
        seqtab=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "dada2",
                "seqTab.nochimeras.RDS"
        ), 
        taxa_table=os.path.join(
                OUTPUT_FOLDER,
                "results",
                "taxa.RDS"
        )
    output:
        count_table=os.path.join(
                OUTPUT_FOLDER,
                "results",
                "count_table_from_dada2.csv"
        ),
        genus_table=os.path.join(
                OUTPUT_FOLDER,
                "results",
                "genus_table_from_dada2.csv"
        ),
        otu_table=os.path.join(
                OUTPUT_FOLDER,
                "results",
                "PICRUSt2_otu_table_with_taxonomy.csv"
        ),
        otu_table_no_taxonomy=os.path.join(
                OUTPUT_FOLDER,
                "results",
                "PICRUSt2_otu_table.csv"
        ),
        repseqs=os.path.join(
                OUTPUT_FOLDER,
                "results",
                "PICRUSt2_representative.csv"
        ),
    log:
        os.path.join(
                OUTPUT_FOLDER,
                "logs",
                "dada2",
                "create_table.log"
        )
    threads: 5 # set desired number of threads here
    conda: 
        "../envs/R.yaml"
    script:
        "../scripts/make_table_dada2.R"

##########################################################################
##########################################################################
