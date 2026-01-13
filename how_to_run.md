# How to Run the DNA Sequence Analysis Workflow

This document provides instructions on how to set up and run the Snakemake workflow for multiple sequence alignment and conservation analysis of DNA sequences.

## 1. Project Structure

The project is organized as follows:

```
.
├── data/
│   └── (Your FASTA files go here)
├── results/
│   └── (Output files will be generated here)
├── config.yaml
├── how_to_run.md
└── Snakefile
```

## 2. Prerequisites

Before running the workflow, ensure you have the following software installed:

*   **Snakemake:** A workflow management system. For this project, Snakemake has been installed into a local Python virtual environment located at `.venv/`.

*   **Docker:** A container platform that allows you to run applications in isolated environments. This workflow uses Docker to ensure reproducibility and easy management of bioinformatics tools (MAFFT and trimAl). Please refer to the [Docker documentation](https://docs.docker.com/get-docker/) for installation instructions. Make sure the Docker daemon is running before executing the workflow.

## 3. Prepare Your Data

1.  **Place FASTA Files:** Put all your DNA sequence files into the `data/` directory.
    *   Each file should be in **FASTA format**.
    *   Ensure all files have a `.fasta` extension (e.g., `isolate1.fasta`, `region_A.fasta`).

## 4. Run the Workflow

1.  **Navigate to Project Directory:** Open your terminal and navigate to the root directory of this project:
    ```bash
    cd /Users/tungchen/Projects/Lmac-DNA
    ```

2.  **Execute Snakemake:** Run the workflow using the following command:

    ```bash
    .venv/bin/snakemake --cores <NUMBER_OF_CORES>
    ```

    *   Replace `<NUMBER_OF_CORES>` with the number of CPU cores you want to allocate to the workflow (e.g., `8`). The `mafft_align` step can utilize multiple threads.
    *   We are running `snakemake` from the `.venv/bin` directory to use the version we installed locally.
    *   The `Snakefile` is now written to automatically use Docker via the `containerized` directive, so the `--use-docker` flag is no longer necessary. Snakemake will automatically pull the required Docker images.

    The workflow will proceed through the following steps:
    *   **Concatenate Sequences:** All `.fasta` files from the `data/` directory will be combined into a single `results/combined.fasta` file.
    *   **Multiple Sequence Alignment (MAFFT):** MAFFT will perform a multiple sequence alignment on the `combined.fasta` file, producing `results/mafft/combined.aligned.fasta`.
    *   **Conservation Analysis (trimAl):** trimAl will then process the aligned FASTA file to remove poorly aligned regions, which can be interpreted as less conserved areas, generating `results/trimal/combined.trimmed.fasta`.

## 5. View Results

The final output, which is the trimmed multiple sequence alignment, will be located at:

*   `results/trimal/combined.trimmed.fasta`

You will also find log files for MAFFT and trimAl in the `results/logs/` directory, which can be useful for debugging or understanding the execution details.
