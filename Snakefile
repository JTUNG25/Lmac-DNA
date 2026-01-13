# Snakefile for Multiple Sequence Alignment and Conservation Analysis

import glob
import os

# Load the configuration
configfile: "config.yaml"

# Get a list of all sample files from the data directory
SAMPLES = glob.glob(os.path.join(config["data_dir"], "*.fasta"))

# The 'all' rule is the default target.
# It specifies the final files we want to create.
rule all:
    input:
        "results/trimal/combined.trimmed.fasta",

# Rule to concatenate all input FASTA files into one
rule concatenate_sequences:
    input:
        SAMPLES
    output:
        "results/combined.fasta"
    shell:
        "cat {input} > {output}"

# Rule for Multiple Sequence Alignment using MAFFT
rule mafft_align:
    input:
        "results/combined.fasta"
    output:
        "results/mafft/combined.aligned.fasta"
    log:
        "results/logs/mafft.log"
    containerized:
        "docker://quay.io/biocontainers/mafft:7.520--h57928b3_1"
    threads: 8
    shell:
        "mafft --thread {threads} --auto {input} > {output} 2> {log}"

# Rule for Conservation Analysis and Trimming with trimAl
rule trimal_analysis:
    input:
        "results/mafft/combined.aligned.fasta"
    output:
        "results/trimal/combined.trimmed.fasta"
    log:
        "results/logs/trimal.log"
    containerized:
        "docker://quay.io/biocontainers/trimal:1.4.1--h57928b3_6"
    shell:
        "trimal -in {input} -out {output} -automated1 > {log}"
