#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --time=04:00:00
#SBATCH --job-name=genomad_viral
#SBATCH --output=genomad.log

# ═════════════════════════════════════════════════════════════════════════════
# geNomad submission script for Bunya HPC
# For viral detection in Trinity assembly
# ═════════════════════════════════════════════════════════════════════════════

# Load conda
source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
mamba activate genomad_env

# Optional: Set temporary directories to avoid NFS bottleneck
export TMPDIR=$TMPDIR

# ─────────────────────────────────────────────────────────────────────────────
# INPUT CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────
TRINITY_FILE="Lmac_Trinity.fasta"          # Your Trinity assembly file
OUTPUT_DIR="viral_results"                 # Output directory name
DATABASE_DIR="genomad_db"            # geNomad database location

# ─────────────────────────────────────────────────────────────────────────────
# RUN geNomad
# ─────────────────────────────────────────────────────────────────────────────

echo "Starting geNomad analysis..."
echo "Input:    $TRINITY_FILE"
echo "Output:   $OUTPUT_DIR"
echo "Database: $DATABASE_DIR"
echo "CPUs:     $SLURM_CPUS_PER_TASK"
echo "Memory:   $SLURM_MEM_PER_NODE MB"
echo "Time:     $(date)"
echo ""

# Run geNomad with conservative filtering (recommended for high confidence)
genomad end-to-end \
    -t $SLURM_CPUS_PER_TASK \
    --conservative \
    "$TRINITY_FILE" \
    "$OUTPUT_DIR" \
    "$DATABASE_DIR"

echo ""
echo "geNomad analysis completed!"
echo "Results saved in: $OUTPUT_DIR/"
echo "Time: $(date)"