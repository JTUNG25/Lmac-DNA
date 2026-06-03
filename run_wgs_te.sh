#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --job-name=wgs_te
#SBATCH --output=wgs_te_sm.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

cd /QRISdata/Q9141/lmac_dna

export TMPDIR=/QRISdata/Q9141/lmac_dna/tmp
export APPTAINER_TMPDIR=$TMPDIR/apptainer_tmp
export APPTAINER_CACHEDIR=$TMPDIR/apptainer_cache
mkdir -p $TMPDIR $APPTAINER_TMPDIR $APPTAINER_CACHEDIR

# ── symlink shared files from RNA project (run once) ─────────────────────────
mkdir -p data/genome data/repeatmodeler

[ ! -f data/genome/JN3.fasta ] && \
    ln -s /QRISdata/Q9141/lmac_rna/data/genome/JN3.fasta \
          data/genome/JN3.fasta

[ ! -f data/genome/JN3.fasta.fai ] && \
    ln -s /QRISdata/Q9141/lmac_rna/data/genome/JN3.fasta.fai \
          data/genome/JN3.fasta.fai

[ ! -f data/genome/JN3.te.gtf ] && \
    ln -s /QRISdata/Q9141/lmac_rna/data/genome/JN3.te.gtf \
          data/genome/JN3.te.gtf

[ ! -f data/genome/JN3.te.bed ] && \
    ln -s /QRISdata/Q9141/lmac_rna/data/genome/JN3.te.bed \
          data/genome/JN3.te.bed

[ ! -f data/repeatmodeler/JN3-families.fa ] && \
    ln -s /QRISdata/Q9141/lmac_rna/results/repeatmodeler/JN3-families.fa \
          data/repeatmodeler/JN3-families.fa

# ── bunya profile (copy from RNA project if not present) ─────────────────────
[ ! -d profiles ] && \
    cp -r /QRISdata/Q9141/lmac_rna/profiles .

chmod +x profiles/bunya/status-sacct-robust.sh

# ── run ───────────────────────────────────────────────────────────────────────
snakemake -s wgs_te.smk --unlock --profile profiles/bunya/

snakemake -s wgs_te.smk \
    --profile profiles/bunya/ \
    --singularity-args "--bind /QRISdata/Q9141 --bind $TMPDIR"