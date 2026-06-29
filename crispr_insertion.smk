#!/usr/bin/env python3
# =============================================================================
# Logic:
#   1. Align WGS reads from CRISPR-mediated mutants to reference genome
#   2. Identify evidence of insertions by detecting:
#      - Discordant read pairs (improperly paired)
#      - Soft-clipped reads (sequences hanging off alignment ends)
#      - Unmapped reads (mate inserted, mate unmapped)
#   3. Extract these reads to separate BAM files for inspection
#   4. Generate hotspot report: positions with clustering of insertion evidence
#   5. QC summary with mapping statistics
# =============================================================================

import os
import re
from pathlib import Path
from collections import defaultdict

bowtie2  = "docker://quay.io/biocontainers/bowtie2:2.5.4--py312h2b63842_1"
samtools = "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
bedtools = "docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"

FASTP_DIR = "results/fastp"  # From wgs_te.smk output
GENOME_DIR = "data/genome"

# ── sample to reference genome mapping ──────────────────────────────────────────

SAMPLE_TO_REFERENCE = {
    # Δago1 group
    "A1-1": f"{GENOME_DIR}/JN3_A1.fa",
    "A1-2": f"{GENOME_DIR}/JN3_A1.fa",
    "A1-3": f"{GENOME_DIR}/JN3_A1.fa",
    # Δago2 group
    "A2-1": f"{GENOME_DIR}/JN3_A2.fa",
    "A2-2": f"{GENOME_DIR}/JN3_A2.fa",
    "A2-3": f"{GENOME_DIR}/JN3_A2.fa",
    # Δago3 group
    "A3-1": f"{GENOME_DIR}/JN3_A3.fa",
    # Δago13 group
    "A13-1": f"{GENOME_DIR}/JN3_A13.fa",
    "A13-2": f"{GENOME_DIR}/JN3_A13.fa",
    # Δdcl1 group
    "D1": f"{GENOME_DIR}/JN3_D1.fa",
    # Δdcl2 group
    "D2-2": f"{GENOME_DIR}/JN3_D2.fa",
    "D2-3": f"{GENOME_DIR}/JN3_D2.fa",
    # Δrdrp1 group
    "R1-1": f"{GENOME_DIR}/JN3_R1.fa",
    "R1-2": f"{GENOME_DIR}/JN3_R1.fa",
    # Δrdrp2 group
    "R2-2": f"{GENOME_DIR}/JN3_R2.fa",
    "R2-3": f"{GENOME_DIR}/JN3_R2.fa",
    "R2-4": f"{GENOME_DIR}/JN3_R2.fa",
    "R2-5": f"{GENOME_DIR}/JN3_R2.fa",
    # Δrdrp3 group
    "R3-1": f"{GENOME_DIR}/JN3_R3.fa",
    "R3-2": f"{GENOME_DIR}/JN3_R3.fa",
    "R3-3": f"{GENOME_DIR}/JN3_R3.fa",
    # Δrdrp12 group
    "R12-1": f"{GENOME_DIR}/JN3_R12.fa",
    "R12-2": f"{GENOME_DIR}/JN3_R12.fa",
    "R12-3": f"{GENOME_DIR}/JN3_R12.fa",
    # WT
    "D5": f"{GENOME_DIR}/JN3.fa",
}

ALL_REFERENCES = sorted(set(SAMPLE_TO_REFERENCE.values()))
ALL_SAMPLES = sorted(SAMPLE_TO_REFERENCE.keys())

# ── helper functions ──────────────────────────────────────────────────────────
def r1(wc): 
    return f"{FASTP_DIR}/{wc.sample}_R1.fastq.gz"

def r2(wc): 
    return f"{FASTP_DIR}/{wc.sample}_R2.fastq.gz"

def get_sample_reference(wc):
    """Return the reference genome for a given sample"""
    return SAMPLE_TO_REFERENCE[wc.sample]


# =============================================================================
# rule all — Final output targets
# =============================================================================
rule all:
    input:
        expand(multiext("{genome}", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
               genome=ALL_REFERENCES),
        # aligned BAMs
        expand("results/bowtie2/{sample}.sorted.bam", sample=ALL_SAMPLES),
        expand("results/bowtie2/{sample}.sorted.bam.bai", sample=ALL_SAMPLES),
        # QC
        expand("results/qc/{sample}.flagstat.txt", sample=ALL_SAMPLES),
        "results/qc/mapping_summary.txt",
        # extracted reads (evidence of insertions)
        expand("results/reads/{sample}.discordant.bam", sample=ALL_SAMPLES),
        expand("results/reads/{sample}.softclipped.bam", sample=ALL_SAMPLES),
        expand("results/reads/{sample}.unmapped_mate.bam", sample=ALL_SAMPLES),
        # insertion hotspot analysis
        expand("results/hotspots/{sample}.insertion_hotspots.bed", sample=ALL_SAMPLES),
        # summary report
        expand("results/summary/{sample}.insertion_candidates.tsv", sample=ALL_SAMPLES),
        "results/summary/all_samples_candidates.tsv",

rule bowtie2_index:

    input:  "{genome}"
    output: 
        multiext("{genome}", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    container: bowtie2
    threads: 8
    resources:
        mem_mb  = 32000,
        runtime = 60,
    params:
        index_name = "{genome}",
    shell:
        "bowtie2-build --threads {threads} {input} {params.index_name}"


rule bowtie2_align:
    """
    Key parameters:
      -a / -k 5  : report multiple alignments (useful for CRISPR repeats)
      -X 1000    : max fragment length (typical insert size range)
      --very-sensitive : thorough but slower
      --no-unal  : exclude unmapped reads from SAM (we extract them separately)
    """
    input:
        r1    = r1,
        r2    = r2,
        ref   = get_sample_reference,
        index = lambda wc: multiext(get_sample_reference(wc), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    output:
        bam = "results/bowtie2/{sample}.sorted.bam",
        bai = "results/bowtie2/{sample}.sorted.bam.bai",
    threads: 16
    resources:
        mem_mb  = 64000,
        runtime = 240,
    params:
        index_name = get_sample_reference,
        rg      = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA",
        tmp     = "/scratch/temp/$SLURM_JOB_ID/bowtie2_{sample}",
        bt2_sif = bowtie2,
        sam_sif = samtools,
    shell:
        """
        mkdir -p {params.tmp}

        apptainer exec {params.bt2_sif} \
            bowtie2 \
                -x {params.index_name} \
                -1 {input.r1} \
                -2 {input.r2} \
                --threads {threads} \
                -k 5 \
                -X 1000 \
                --very-sensitive \
                --no-unal \
                -R '{params.rg}' | \
        apptainer exec {params.sam_sif} \
            samtools sort \
                -@ {threads} \
                -T {params.tmp}/sort \
                -o {output.bam}

        apptainer exec {params.sam_sif} \
            samtools index -@ {threads} {output.bam}

        rm -rf {params.tmp}
        """


# =============================================================================
# SECTION 2 — QC AND MAPPING STATISTICS
# =============================================================================

rule flagstat:
    """
    Generate mapping statistics for each sample.
    Shows: total reads, mapped reads, paired reads, properly paired, etc.
    """
    input:
        bam = "results/bowtie2/{sample}.sorted.bam",
        bai = "results/bowtie2/{sample}.sorted.bam.bai",
    output: 
        "results/qc/{sample}.flagstat.txt"
    container: samtools
    threads: 2
    resources:
        mem_mb  = 8000,
        runtime = 15,
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output}"


rule mapping_summary:
    """
    Summarize mapping quality across all samples.
    Check for unexpectedly low mapped reads or high multi-mapping rates.
    """
    input:
        expand("results/qc/{sample}.flagstat.txt", sample=ALL_SAMPLES),
    output:
        "results/qc/mapping_summary.txt",
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 5,
    run:
        with open(output[0], "w") as out:
            out.write("sample\ttotal_reads\tmapped_reads\tmapped_pct\t"
                     "properly_paired\tproperly_paired_pct\tnote\n")
            for sample in sorted(ALL_SAMPLES):
                with open(f"results/qc/{sample}.flagstat.txt") as f:
                    lines = f.readlines()
                
                total    = int(lines[0].split()[0])
                mapped   = int([l for l in lines if "mapped (" in l][0].split()[0])
                properly = int([l for l in lines if "properly paired" in l][0].split()[0])
                
                mapped_pct    = f"{100*mapped/total:.1f}" if total > 0 else "0"
                properly_pct  = f"{100*properly/total:.1f}" if total > 0 else "0"
                
                note = ""
                if total < 1_000_000:
                    note = "WARNING: low read count"
                elif mapped/total < 0.7:
                    note = "WARNING: low mapping rate"
                
                out.write(f"{sample}\t{total}\t{mapped}\t{mapped_pct}%\t"
                         f"{properly}\t{properly_pct}%\t{note}\n")


# =============================================================================
# SECTION 3 — EXTRACT READS WITH INSERTION EVIDENCE
#
# FLAG meanings:
#   2  = properly paired (both mates aligned as expected)
#   4  = unmapped
#   8  = mate unmapped
#   16 = reverse strand
#   32 = mate on reverse strand
#
# Discordant (use -F 2):
#   - Pairs far apart
#   - Unexpected orientation (both same strand)
#   - Different chromosomes
#   All suggest structural variation / insertion
#
# Soft-clipped (CIGAR has S):
#   - Reads with unaligned sequence at ends
#   - Common at insertion junctions
#
# Unmapped with mapped mate (use -f 4 -F 8):
#   - Insert broke one side of pair
#   - Other mate still aligns nearby
# =============================================================================

rule extract_discordant:
    """
    Extract improperly paired reads (discordant pairs).
    These are pairs where:
    - One or both mates are unmapped
    - Mates are on different chromosomes
    - Mates are far apart (>1000bp default)
    - Mates are in unexpected orientation
    
    All of these suggest structural variation including insertions.
    
    Flag -F 2 = exclude flag 2 (properly paired) → keep only improper pairs
    """
    input:
        bam = "results/bowtie2/{sample}.sorted.bam",
        bai = "results/bowtie2/{sample}.sorted.bam.bai",
    output:
        "results/reads/{sample}.discordant.bam",
    container: samtools
    threads: 4
    resources:
        mem_mb  = 16000,
        runtime = 30,
    shell:
        """
        samtools view \
            -b \
            -@ {threads} \
            -F 2 \
            {input.bam} > {output}
        """


rule extract_softclipped:
    """
    Extract soft-clipped reads.
    
    Soft clipping (CIGAR S) = part of read is aligned, part is not.
    Common at insertion breakpoints where the junction sequence 
    doesn't match the reference.
    
    Extract reads with any soft-clipping in CIGAR string.
    """
    input:
        bam = "results/bowtie2/{sample}.sorted.bam",
        bai = "results/bowtie2/{sample}.sorted.bam.bai",
    output:
        "results/reads/{sample}.softclipped.bam",
    container: samtools
    threads: 4
    resources:
        mem_mb  = 16000,
        runtime = 30,
    run:
        import pysam
        
        # Open BAM and create output BAM with same header
        inbam = pysam.AlignmentFile(input.bam, "rb")
        outbam = pysam.AlignmentFile(output[0], "wb", template=inbam)
        
        count = 0
        for read in inbam:
            # Check if CIGAR contains S (soft-clipping)
            if read.cigarstring and "S" in read.cigarstring:
                outbam.write(read)
                count += 1
        
        inbam.close()
        outbam.close()
        
        print(f"  Extracted {count} soft-clipped reads from {input.bam}")


rule extract_unmapped_mate:
    """
    Extract reads where one mate is unmapped but the other is mapped.
    
    Flag -f 4  = include unmapped reads
    Flag -F 8  = exclude reads where mate is also unmapped
    
    Result: reads that are unmapped but have a mapped mate
    (suggests insert at this location broke one side of the pair)
    """
    input:
        bam = "results/bowtie2/{sample}.sorted.bam",
        bai = "results/bowtie2/{sample}.sorted.bam.bai",
    output:
        "results/reads/{sample}.unmapped_mate.bam",
    container: samtools
    threads: 4
    resources:
        mem_mb  = 16000,
        runtime = 30,
    shell:
        """
        samtools view \
            -b \
            -@ {threads} \
            -f 4 \
            -F 8 \
            {input.bam} > {output}
        """


# =============================================================================
# SECTION 4 — IDENTIFY INSERTION HOTSPOTS
#
# For each read type (discordant, soft-clipped, unmapped_mate),
# extract the mapping positions and count coverage in sliding windows.
# Regions with clustering of insertion evidence are candidates for
# actual insertion sites.
#
# Output: BED file with positions and evidence counts
# =============================================================================

rule find_hotspots:
    """
    Find insertion hotspots by analyzing clustering of evidence reads.
    
    For each sample, combines:
    - Positions where discordant pairs map
    - Positions of soft-clipped reads
    - Positions of unmapped reads (mate)
    
    Produces a BED file with hotspot regions and evidence counts.
    """
    input:
        discordant = "results/reads/{sample}.discordant.bam",
        softclipped = "results/reads/{sample}.softclipped.bam",
        unmapped   = "results/reads/{sample}.unmapped_mate.bam",
        fai        = lambda wc: get_sample_reference(wc) + ".fai",
    output:
        hotspots = "results/hotspots/{sample}.insertion_hotspots.bed",
    container: samtools
    threads: 2
    resources:
        mem_mb  = 16000,
        runtime = 30,
    run:
        import pysam
        from collections import defaultdict
        
        # Collect evidence positions
        positions = defaultdict(lambda: {"discordant": 0, "softclipped": 0, "unmapped": 0})
        
        # Read discordant pairs
        try:
            bam = pysam.AlignmentFile(input.discordant, "rb")
            for read in bam:
                if not read.is_unmapped:
                    key = f"{read.reference_name}:{read.reference_start}"
                    positions[key]["discordant"] += 1
            bam.close()
        except:
            pass
        
        # Read soft-clipped
        try:
            bam = pysam.AlignmentFile(input.softclipped, "rb")
            for read in bam:
                if not read.is_unmapped:
                    key = f"{read.reference_name}:{read.reference_start}"
                    positions[key]["softclipped"] += 1
            bam.close()
        except:
            pass
        
        # Read unmapped with mate
        try:
            bam = pysam.AlignmentFile(input.unmapped, "rb")
            for read in bam:
                # Get mate position if available
                if read.next_reference_name and not read.next_reference_name.startswith("*"):
                    key = f"{read.next_reference_name}:{read.next_reference_start}"
                    positions[key]["unmapped"] += 1
            bam.close()
        except:
            pass
        
        # Write hotspots (BED format)
        with open(output.hotspots, "w") as out:
            out.write("# Insertion evidence hotspots\n")
            out.write("# chrom\tstart\tend\tevidence_score\tdiscordant_count\t"
                     "softclipped_count\tunmapped_mate_count\n")
            
            for key in sorted(positions):
                chrom, pos = key.split(":")
                pos = int(pos)
                evidence = positions[key]
                
                # Simple scoring: sum of all evidence types
                score = (evidence["discordant"] * 2 +  # weight discordant more
                        evidence["softclipped"] * 1 +
                        evidence["unmapped"] * 1)
                
                # Only report if there's meaningful evidence
                if score >= 2:
                    out.write(f"{chrom}\t{pos}\t{pos+1}\t{score}\t"
                             f"{evidence['discordant']}\t"
                             f"{evidence['softclipped']}\t"
                             f"{evidence['unmapped']}\n")


# =============================================================================
# SECTION 5 — GENERATE INSERTION CANDIDATE SUMMARY
#
# For each sample, create a table of:
# - Hotspot location (chrom:pos)
# - Number of supporting reads by type
# - Coverage at that position
# - Ranked by evidence strength
#
# Ready for IGV inspection or PCR validation
# =============================================================================

rule candidate_summary:
    """
    Create a summary table of insertion candidates per sample.
    Ranks hotspots by evidence strength and includes genome context.
    """
    input:
        hotspots = "results/hotspots/{sample}.insertion_hotspots.bed",
        bam      = "results/bowtie2/{sample}.sorted.bam",
        bai      = "results/bowtie2/{sample}.sorted.bam.bai",
    output:
        "results/summary/{sample}.insertion_candidates.tsv",
    container: samtools
    threads: 2
    resources:
        mem_mb  = 8000,
        runtime = 15,
    run:
        import pysam
        
        # Read hotspots
        hotspots_list = []
        with open(input.hotspots) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 7:
                    chrom, start, end, score, disc, soft, unmapped = parts[:7]
                    hotspots_list.append({
                        "chrom": chrom,
                        "pos": int(start),
                        "score": int(score),
                        "discordant": int(disc),
                        "softclipped": int(soft),
                        "unmapped": int(unmapped),
                    })
        
        # Sort by evidence score
        hotspots_list.sort(key=lambda x: x["score"], reverse=True)
        
        # Get coverage at each position
        bam = pysam.AlignmentFile(input.bam, "rb")
        
        with open(output[0], "w") as out:
            out.write("rank\tchrom\tposition\ttotal_evidence_score\t"
                     "discordant_reads\tsoftclipped_reads\tunmapped_mate_reads\t"
                     "local_coverage\tpriority\n")
            
            for rank, hs in enumerate(hotspots_list, 1):
                # Get coverage at position
                try:
                    coverage = bam.count(hs["chrom"], hs["pos"], hs["pos"] + 1)
                except:
                    coverage = "NA"
                
                # Assign priority based on evidence
                if hs["score"] >= 10 and hs["discordant"] >= 3:
                    priority = "HIGH"
                elif hs["score"] >= 5 and hs["softclipped"] >= 2:
                    priority = "MEDIUM"
                else:
                    priority = "LOW"
                
                out.write(f"{rank}\t{hs['chrom']}\t{hs['pos']}\t{hs['score']}\t"
                         f"{hs['discordant']}\t{hs['softclipped']}\t"
                         f"{hs['unmapped']}\t{coverage}\t{priority}\n")
        
        bam.close()


rule merge_all_candidates:
    """
    Merge insertion candidates from all samples into one master table.
    Adds sample ID for cross-sample comparison.
    """
    input:
        expand("results/summary/{sample}.insertion_candidates.tsv", 
               sample=ALL_SAMPLES),
    output:
        "results/summary/all_samples_candidates.tsv",
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 5,
    run:
        header_written = False
        with open(output[0], "w") as out:
            out.write("sample\t")  # Add sample column
            for fp in sorted(input):
                sample = fp.split("/")[-1].replace(".insertion_candidates.tsv", "")
                with open(fp) as f:
                    for line in f:
                        if line.startswith("rank"):
                            if not header_written:
                                out.write(line)
                                header_written = True
                        elif not line.startswith("#"):
                            out.write(f"{sample}\t{line}")
