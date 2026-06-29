#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --job-name=bowtie2
#SBATCH --output=bt.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

export TMPDIR=/QRISdata/Q9141/lmac_dna/tmp
mkdir -p $TMPDIR

# Make status script executable
chmod +x profiles/bunya/status-sacct-robust.sh

# Run snakemake
snakemake -s crispr_insertion.smk --profile profiles/bunya/ 
