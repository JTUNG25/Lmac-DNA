# Snakefile for Multiple Sequence Alignment and Conservation Analysis
# Works on both Mac (with conda) and Linux (with singularity/containers)

import glob
import os

# Load the configuration
configfile: "config_test.yaml"

# Detect platform (Mac vs Linux)
PLATFORM = "mac" if os.uname().sysname == "Darwin" else "linux"

# Get a list of sample files
if config.get("test_mode"):
    test_files = config["test_samples"]
    SAMPLES = [os.path.join(config["data_dir"], s) for s in test_files]
else:
    SAMPLES = glob.glob(os.path.join(config["data_dir"], "*.tar.gz"), recursive=True)

# Extract sample names without extension for use in rules
SAMPLE_NAMES = [os.path.basename(s).replace(".tar.gz", "") for s in SAMPLES]

rule all:
    input:
        "results/trimal/combined.trimmed.fasta",
        "results/snp-sites/polymorphic_sites.vcf"

# Rule to extract tar.gz files
rule extract_samples:
    input:
        lambda wildcards: [s for s in SAMPLES if os.path.basename(s).replace(".tar.gz", "") == wildcards.sample][0]
    output:
        "results/extracted/{sample}.fasta"
    shell:
        "tar -xzf {input} -O > {output}"

# Rule to concatenate all extracted FASTA files into one
rule concatenate_sequences:
    input:
        expand("results/extracted/{sample}.fasta", sample=SAMPLE_NAMES)
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
    conda:
        "environment.yaml" if PLATFORM == "mac" else None
    container:
        "docker://quay.io/biocontainers/mafft:7.520--h57928b3_1" if PLATFORM == "linux" else None
    threads: 2
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
    conda:
        "environment.yaml" if PLATFORM == "mac" else None
    container:
        "docker://quay.io/biocontainers/trimal:1.4.1--h57928b3_6" if PLATFORM == "linux" else None
    shell:
        "trimal -in {input} -out {output} -automated1 > {log}"

# Rule for identifying polymorphic sites using SNP-sites
rule snp_sites_analysis:
    input:
        "results/mafft/combined.aligned.fasta"
    output:
        "results/snp-sites/polymorphic_sites.vcf"
    log:
        "results/logs/snp-sites.log"
    conda:
        "environment.yaml" if PLATFORM == "mac" else None
    container:
        "docker://quay.io/biocontainers/snp-sites:2.5.1--h103a89f_4" if PLATFORM == "linux" else None
    shell:
        "snp-sites -v -o {output} {input} 2> {log}"