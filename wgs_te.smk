#!/usr/bin/env python3
# =============================================================================
# wgs_te.smk  —  WGS TE copy-number screen
# Leptosphaeria maculans JN3
#
# Logic:
#   1. Read significant_TEs.csv (from R script) → candidate TE family list
#      Only upregulated TEs: log2FC > 1, padj < 0.05, already filtered for
#      simple/low-complexity repeats by R script
#   2. Filter JN3.te.bed to candidate families only
#   3. fastp  → trim WGS reads
#   4. BWA-MEM → align all samples to plain JN3.fasta (same reference for all)
#   5. bedtools coverage → mean read depth per TE locus
#   6. samtools coverage → genome-wide mean depth (for normalisation)
#   7. Collapse loci → per-family normalised depth per sample
#   8. Compare each mutant vs WT → log2FC, flag INCREASED if >= log2(1.5)
#   9. Final table cross-referenced against RNA log2FC from significant_TEs.csv
#
# WT: D5 (single replicate — depth comparisons are exploratory screen only)
# No statistical test on depth — biological replicates needed for that.
# =============================================================================

import os
import math
import collections
from pathlib import Path

# ── containers ────────────────────────────────────────────────────────────────
fastp    = "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_3"
bwa      = "docker://quay.io/biocontainers/bwa:0.7.18--he4a0461_0"
samtools = "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
bedtools = "docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"

# ── paths ─────────────────────────────────────────────────────────────────────
RAW_DIR   = "/QRISdata/Q9141/lmac_dna/merged_data"
GENOME    = "data/genome/JN3.fasta"
TE_BED    = "data/genome/JN3.te.bed"
TE_GTF    = "data/genome/JN3.te.gtf"
SIG_TE    = "data/significant_TEs.csv"   # from R script, already on HPC

# ── samples ───────────────────────────────────────────────────────────────────
WGS_WT      = "D5"
WGS_MUTANTS = sorted([
    p.name.replace("_R1.fastq.gz", "")
    for p in Path(RAW_DIR).glob("*_R1.fastq.gz")
    if p.name.replace("_R1.fastq.gz", "") != WGS_WT
])
ALL_WGS = WGS_MUTANTS + [WGS_WT]

# ── mutant → R-script label mapping ──────────────────────────────────────────
# Links WGS sample names to the mutant labels in significant_TEs.csv
# so we can pull the right candidate TE list per mutant group
WGS_TO_LABEL = {
    "A1-1":  "Δago1",  "A1-2":  "Δago1",  "A1-3":  "Δago1",
    "A3-1":  "Δago3",
    "A13-1": "Δago13", "A13-2": "Δago13",
    "D1":    "Δdcl1",
    "D2-2":  "Δdcl2",  "D2-3":  "Δdcl2",
    "R1-1":  "Δrdrp1",
    "R2-2":  "Δrdrp2", "R2-3":  "Δrdrp2",
    "R2-4":  "Δrdrp2", "R2-5":  "Δrdrp2",
    "R12-1": "Δrdrp12","R12-2": "Δrdrp12","R12-3": "Δrdrp12",
}

# ── mutant groups (for per-group depth comparison output) ─────────────────────
GROUPS = {
    "ago1":  {"label": "Δago1",   "samples": ["A1-1", "A1-2", "A1-3"]},
    "ago3":  {"label": "Δago3",   "samples": ["A3-1"]},
    "ago13": {"label": "Δago13",  "samples": ["A13-1", "A13-2"]},
    "dcl1":  {"label": "Δdcl1",   "samples": ["D1"]},
    "dcl2":  {"label": "Δdcl2",   "samples": ["D2-2", "D2-3"]},
    "rdrp1": {"label": "Δrdrp1",  "samples": ["R1-1"]},
    "rdrp2": {"label": "Δrdrp2",  "samples": ["R2-2", "R2-3", "R2-4", "R2-5"]},
    "rdrp12":{"label": "Δrdrp12", "samples": ["R12-1", "R12-2", "R12-3"]},
}

# ── helper ────────────────────────────────────────────────────────────────────
def r1(wc): return f"{RAW_DIR}/{wc.sample}_R1.fastq.gz"
def r2(wc): return f"{RAW_DIR}/{wc.sample}_R2.fastq.gz"


# =============================================================================
# rule all
# =============================================================================
rule all:
    input:
        # trimmed reads
        expand("results/fastp/{sample}_R1.fastq.gz",  sample=ALL_WGS),
        # aligned BAMs
        expand("results/bwa/{sample}.sorted.bam",     sample=ALL_WGS),
        expand("results/bwa/{sample}.sorted.bam.bai", sample=ALL_WGS),
        # QC
        expand("results/qc/{sample}.flagstat.txt",    sample=ALL_WGS),
        "results/qc/coverage_summary.txt",
        # candidate BED (TE loci for upregulated families only)
        "results/te_depth/candidate_te_loci.bed",
        # per-sample depth
        expand("results/te_depth/{sample}.locus_depth.bed",   sample=ALL_WGS),
        expand("results/te_depth/{sample}.genome_depth.txt",  sample=ALL_WGS),
        expand("results/te_depth/{sample}.family_depth.tsv",  sample=ALL_WGS),
        # per-group comparison vs WT
        expand("results/te_depth/{group}_vs_wt.tsv", group=GROUPS),
        # final cross-referenced summary
        "results/crossref/all_groups_summary.tsv",


# =============================================================================
# SECTION 0 — FASTP TRIMMING
# =============================================================================

rule fastp_trim:
    input:
        r1 = r1,
        r2 = r2,
    output:
        r1   = "results/fastp/{sample}_R1.fastq.gz",
        r2   = "results/fastp/{sample}_R2.fastq.gz",
        html = "results/fastp/{sample}_fastp.html",
        json = "results/fastp/{sample}_fastp.json",
    container: fastp
    threads: 8
    resources:
        mem_mb  = 16000,
        runtime = 120,
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --html {output.html} --json {output.json} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --correction \
            --qualified_quality_phred 20 \
            --length_required 50
        """


# =============================================================================
# SECTION 1 — BWA-MEM ALIGNMENT
# All 18 samples aligned to the same plain JN3.fasta.
# Critical: no T-DNA reference, so depth is comparable across all samples.
# BWA-MEM chosen over Bowtie2 because -a flag retains all alignments for
# multi-mapping reads, which is important for repetitive TE sequences.
# Piped directly into samtools sort — no intermediate SAM file on disk.
# =============================================================================

rule bwa_index:
    input:  GENOME
    output: multiext(GENOME, ".amb", ".ann", ".bwt", ".pac", ".sa")
    container: bwa
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 30,
    shell:
        "bwa index {input}"


rule bwa_mem_align:
    input:
        r1    = "results/fastp/{sample}_R1.fastq.gz",
        r2    = "results/fastp/{sample}_R2.fastq.gz",
        index = multiext(GENOME, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam = "results/bwa/{sample}.sorted.bam",
        bai = "results/bwa/{sample}.sorted.bam.bai",
    threads: 16
    resources:
        mem_mb  = 64000,
        runtime = 240,
    params:
        genome  = GENOME,
        rg      = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA",
        tmp     = "/scratch/temp/$SLURM_JOB_ID/bwa_{sample}",
        bwa_sif = bwa,
        sam_sif = samtools,
    shell:
        """
        mkdir -p {params.tmp}

        apptainer exec {params.bwa_sif} \
            bwa mem \
                -t {threads} \
                -a \
                -M \
                -R '{params.rg}' \
                {params.genome} \
                {input.r1} {input.r2} | \
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
# SECTION 2 — QC
# =============================================================================

rule flagstat:
    input:
        bam = "results/bwa/{sample}.sorted.bam",
        bai = "results/bwa/{sample}.sorted.bam.bai",
    output: "results/qc/{sample}.flagstat.txt"
    container: samtools
    threads: 2
    resources:
        mem_mb  = 4000,
        runtime = 15,
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output}"


rule coverage_summary:
    """
    One-line summary per sample: total reads, mapped reads, mapping %.
    Check this before interpreting depth results — any sample with
    unexpectedly low read counts should be investigated.
    """
    input:
        expand("results/qc/{sample}.flagstat.txt", sample=ALL_WGS),
    output:
        "results/qc/coverage_summary.txt",
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = 5,
    run:
        with open(output[0], "w") as out:
            out.write("sample\ttotal_reads\tmapped_reads\tmapping_pct\tnote\n")
            for sample in ALL_WGS:
                with open(f"results/qc/{sample}.flagstat.txt") as f:
                    lines = f.readlines()
                total  = int(lines[0].split()[0])
                mapped = int([l for l in lines if "mapped (" in l][0].split()[0])
                pct    = f"{100*mapped/total:.1f}" if total > 0 else "0"
                note   = "WARNING: low read count — check before interpreting" \
                         if total < 1_000_000 else ""
                out.write(f"{sample}\t{total}\t{mapped}\t{pct}%\t{note}\n")


# =============================================================================
# SECTION 3 — BUILD CANDIDATE TE LOCI BED
#
# Reads significant_TEs.csv (already filtered by R script:
#   - padj < 0.05, log2FC > 1 (upregulated only)
#   - simple/low-complexity repeats removed
#   - gene rows removed)
# Extracts the te_name from the feature column (everything before first ":")
# Filters JN3.te.bed to only loci belonging to those candidate families.
#
# This means depth is only computed over TE families you actually care about,
# not the entire repeat content of the genome.
# =============================================================================

rule make_candidate_bed:
    """
    Extract candidate TE family names from significant_TEs.csv and
    filter JN3.te.bed to those families only.

    BED col4 format: locus_id;te_name  (e.g. JN3_BAC01_1074_1337;rnd-1_family-13)
    We match on the te_name part (after the semicolon).
    """
    input:
        sig_csv = SIG_TE,
        te_bed  = TE_BED,
    output:
        "results/te_depth/candidate_te_loci.bed",
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 5,
    run:
        # collect all upregulated te_names across all mutants
        candidate_names = set()
        with open(input.sig_csv) as f:
            f.readline()  # skip header
            for line in f:
                parts  = line.strip().split(",")
                if len(parts) < 2:
                    continue
                feature = parts[1].strip('"')
                te_name = feature.split(":")[0]   # e.g. rnd-1_family-29
                candidate_names.add(te_name)

        print(f"  {len(candidate_names)} candidate TE families from R script")

        # filter BED to matching loci
        kept = 0
        with open(input.te_bed) as f_in, open(output[0], "w") as f_out:
            for line in f_in:
                if line.startswith("#") or not line.strip():
                    continue
                col4    = line.strip().split("\t")[3]   # locus_id;te_name
                te_name = col4.split(";")[1] if ";" in col4 else col4
                if te_name in candidate_names:
                    f_out.write(line)
                    kept += 1

        print(f"  {kept} TE loci kept in candidate BED")


# =============================================================================
# SECTION 4 — TE READ DEPTH
#
# How copy number is estimated from WGS read depth:
#
#   The reference genome has known TE loci (positions in JN3.te.bed).
#   After BWA-MEM alignment, reads pile up over these loci.
#   More genomic copies of a TE family → more reads mapping to those positions
#   → higher read depth at those loci.
#
#   Raw depth varies between samples due to differences in total sequencing
#   depth, so we normalise:
#
#     norm_depth = mean_depth_over_TE_loci / genome_wide_mean_depth
#
#   This gives a value independent of sequencing depth.
#   Then: log2FC = log2(mutant_norm / wt_norm)
#   A positive log2FC suggests more TE copies in the mutant vs WT.
#
#   Important caveat: this measures depth at REFERENCE loci only.
#   New insertions at novel sites are not captured here — they would
#   require McClintock-style split-read analysis (can be added later).
# =============================================================================

rule te_locus_depth:
    """
    Mean read depth over each candidate TE locus.
    bedtools coverage -mean: for each interval in the BED, compute the
    mean per-base depth from the BAM. Output adds one column (mean depth)
    to the BED.
    Only runs over candidate loci (not the full TE annotation) for efficiency.
    """
    input:
        bam = "results/bwa/{sample}.sorted.bam",
        bai = "results/bwa/{sample}.sorted.bam.bai",
        bed = "results/te_depth/candidate_te_loci.bed",
    output:
        "results/te_depth/{sample}.locus_depth.bed",
    container: bedtools
    threads: 2
    resources:
        mem_mb  = 16000,
        runtime = 60,
    shell:
        """
        bedtools coverage \
            -a {input.bed} \
            -b {input.bam} \
            -mean \
        > {output}
        """


rule genome_mean_depth:
    """
    Genome-wide mean depth for normalisation.
    samtools coverage outputs per-contig stats; we compute a weighted mean
    across all contigs (weighted by contig length) to get one value per sample.
    This is the denominator in the normalisation formula above.
    """
    input:
        bam = "results/bwa/{sample}.sorted.bam",
        bai = "results/bwa/{sample}.sorted.bam.bai",
    output:
        "results/te_depth/{sample}.genome_depth.txt",
    container: samtools
    threads: 2
    resources:
        mem_mb  = 8000,
        runtime = 20,
    shell:
        """
        samtools coverage {input.bam} | \
        awk 'NR>1 && $3>0 {{
            bases += $3
            depth += $7 * $3
        }} END {{
            printf "%.4f\\n", (bases>0 ? depth/bases : 0)
        }}' > {output}
        """


rule family_depth:
    """
    Collapse per-locus depth → per-family mean normalised depth.

    For each TE family, there may be many individual loci (copies) in the
    reference genome. We average the depth across all loci belonging to the
    same family, then divide by genome-wide mean depth.

    This gives one normalised value per family per sample, directly comparable
    to the log2FC values from TEtranscripts — just from the DNA side.
    """
    input:
        depth  = "results/te_depth/{sample}.locus_depth.bed",
        gd     = "results/te_depth/{sample}.genome_depth.txt",
        te_gtf = TE_GTF,
    output:
        "results/te_depth/{sample}.family_depth.tsv",
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 10,
    run:
        import re, collections

        with open(input.gd) as f:
            gd = float(f.read().strip() or 0)

        # family → class from GTF (for annotation in output)
        fam_class = {}
        with open(input.te_gtf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                mg = re.search(r'gene_id "([^"]+)"', line)
                mc = re.search(r'class_id "([^"]+)"', line)
                if mg and mc:
                    fam_class[mg.group(1)] = mc.group(1)

        # accumulate per-locus depths by family
        # BED col4: locus_id;te_name  → split on ";" to get te_name
        fam_depths = collections.defaultdict(list)
        n_loci     = collections.defaultdict(int)
        with open(input.depth) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                col4      = parts[3]
                te_name   = col4.split(";")[1] if ";" in col4 else col4
                locus_dep = float(parts[6])   # mean depth from bedtools -mean
                fam_depths[te_name].append(locus_dep)
                n_loci[te_name] += 1

        with open(output[0], "w") as out:
            out.write(
                "# Normalised TE family depth — EXPLORATORY SCREEN\n"
                "# norm_depth = mean_locus_depth / genome_mean_depth\n"
                "te_name\tclass\tsample\tn_loci\t"
                "mean_raw_depth\tgenome_mean_depth\tnorm_depth\n"
            )
            for te_name in sorted(fam_depths):
                depths = fam_depths[te_name]
                mean_d = sum(depths) / len(depths)
                norm_d = mean_d / gd if gd > 0 else 0
                cls    = fam_class.get(te_name, "Unknown")
                out.write(
                    f"{te_name}\t{cls}\t{wildcards.sample}\t"
                    f"{n_loci[te_name]}\t{mean_d:.4f}\t{gd:.4f}\t{norm_d:.4f}\n"
                )


# =============================================================================
# SECTION 5 — COMPARE MUTANT DEPTH VS WT
#
# For each mutant group:
#   log2FC_depth = log2(mutant_norm_depth / wt_norm_depth)
#
#   INCREASED flag: log2FC_depth >= log2(1.5)  i.e. >= 1.5-fold more reads
#                   over TE loci in mutant vs WT
#   NEW flag:       TE family present in mutant but zero depth in WT
#                   (rare but biologically interesting)
#
# One column per WGS sample so you can see if the signal is consistent
# within a mutant group (e.g. all three Δago1 samples show INCREASED).
# =============================================================================

rule depth_vs_wt:
    input:
        wt      = f"results/te_depth/{WGS_WT}.family_depth.tsv",
        mutants = lambda wc: expand(
            "results/te_depth/{sample}.family_depth.tsv",
            sample=GROUPS[wc.group]["samples"]
        ),
    output:
        "results/te_depth/{group}_vs_wt.tsv",
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 5,
    run:
        import math, re

        group   = wildcards.group
        members = GROUPS[group]["samples"]

        def parse_depth(path):
            result = {}
            with open(path) as f:
                for line in f:
                    if line.startswith("#") or line.startswith("te_name"):
                        continue
                    p = line.strip().split("\t")
                    if len(p) < 7:
                        continue
                    # te_name, class, sample, n_loci, mean_raw, genome_mean, norm
                    result[p[0]] = {
                        "class":  p[1],
                        "n_loci": int(p[3]),
                        "norm":   float(p[6]),
                    }
            return result

        wt_data  = parse_depth(input.wt)
        mut_data = {}
        for fp in input.mutants:
            sname = re.search(r"te_depth/(.+)\.family_depth", fp).group(1)
            mut_data[sname] = parse_depth(fp)

        all_te = set(wt_data)
        for d in mut_data.values():
            all_te |= set(d)

        mut_hdr = "\t".join(
            [f"{s}_norm_depth\t{s}_log2FC\t{s}_flag" for s in members]
        )

        with open(output[0], "w") as out:
            out.write(
                "# EXPLORATORY SCREEN — single WT replicate, no statistical testing\n"
                "# log2FC_depth = log2(mutant_norm / wt_norm)\n"
                "# INCREASED: log2FC >= 0.585 (>=1.5-fold vs WT)\n"
                "# NEW: TE present in mutant, zero depth in WT\n"
                f"te_name\tclass\tn_reference_loci\twt_norm_depth\t{mut_hdr}\n"
            )

            for te in sorted(all_te):
                wt_info  = wt_data.get(te, {"class": "Unknown", "n_loci": 0, "norm": 0.0})
                wt_norm  = wt_info["norm"]
                cls      = wt_info["class"]
                n_loci   = wt_info["n_loci"]
                row = [te, cls, str(n_loci), f"{wt_norm:.4f}"]

                for samp in members:
                    mut_norm = mut_data[samp].get(te, {}).get("norm", 0.0)
                    if wt_norm > 0 and mut_norm > 0:
                        lfc  = math.log2(mut_norm / wt_norm)
                        flag = "INCREASED" if lfc >= math.log2(1.5) else ""
                    elif mut_norm > 0:
                        lfc  = float("inf")
                        flag = "NEW"
                    else:
                        lfc  = float("-inf")
                        flag = ""
                    lfc_str = f"{lfc:.4f}" if math.isfinite(lfc) else str(lfc)
                    row += [f"{mut_norm:.4f}", lfc_str, flag]

                out.write("\t".join(row) + "\n")


# =============================================================================
# SECTION 6 — FINAL CROSS-REFERENCE SUMMARY
#
# Joins per-group depth results with the RNA log2FC from significant_TEs.csv
# so the final table has both RNA and DNA evidence side by side.
#
# Output columns per TE family:
#   te_name, te_family, te_class        — identity
#   rna_log2FC, rna_padj               — from TEtranscripts via R script
#   n_reference_loci                   — how many copies in reference genome
#   wt_norm_depth                      — normalised depth in WT
#   [per sample] norm_depth, log2FC, flag  — DNA depth evidence
#   any_increased                      — TRUE if any sample flagged INCREASED
#   consistent_increased               — TRUE if ALL samples in group INCREASED
#                                        (stronger evidence when >1 sample)
#
# Evidence interpretation:
#   RNA up + depth INCREASED in all samples  → strong candidate
#   RNA up + depth INCREASED in some samples → moderate candidate
#   RNA up + depth flat                      → transcriptional derepression only
#                                              (no copy number change detected)
# =============================================================================

rule crossref_group:
    input:
        sig_csv   = SIG_TE,
        depth_tsv = "results/te_depth/{group}_vs_wt.tsv",
    output:
        "results/crossref/{group}_summary.tsv",
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 5,
    run:
        import math

        group   = wildcards.group
        label   = GROUPS[group]["label"]
        members = GROUPS[group]["samples"]

        # ── RNA sig results for this mutant label ─────────────────────────
        rna_data = {}   # te_name -> (log2FC, padj, te_family, te_class)
        with open(input.sig_csv) as f:
            f.readline()
            for line in f:
                parts   = line.strip().split(",")
                if len(parts) < 6:
                    continue
                mutant  = parts[0].strip('"')
                if mutant != label:
                    continue
                feature = parts[1].strip('"')
                fp      = feature.split(":")
                te_name = fp[0]
                te_fam  = fp[1] if len(fp) > 1 else "Unknown"
                te_cls  = fp[2] if len(fp) > 2 else "Unknown"
                try:
                    lfc  = float(parts[2])
                    padj = float(parts[5])
                except ValueError:
                    continue
                rna_data[te_name] = (lfc, padj, te_fam, te_cls)

        # ── depth results ─────────────────────────────────────────────────
        # parse the group depth file — already has one row per te_name
        depth_rows = {}
        with open(input.depth_tsv) as f:
            header = None
            for line in f:
                if line.startswith("#"):
                    continue
                if line.startswith("te_name"):
                    header = line.strip().split("\t")
                    continue
                parts = line.strip().split("\t")
                if parts:
                    depth_rows[parts[0]] = parts   # keyed by te_name

        # ── write output ──────────────────────────────────────────────────
        samp_hdr = "\t".join(
            [f"{s}_norm_depth\t{s}_log2FC_depth\t{s}_depth_flag"
             for s in members]
        )
        with open(output[0], "w") as out:
            out.write(
                "# Cross-reference: RNA upregulation (TEtranscripts) vs "
                "WGS read depth\n"
                "# EXPLORATORY SCREEN — validate candidates with replicates\n"
            )
            out.write(
                f"group\tmutant\tte_name\tte_family\tte_class\t"
                f"rna_log2FC\trna_padj\t"
                f"n_reference_loci\twt_norm_depth\t"
                f"{samp_hdr}\t"
                f"any_depth_increased\tconsistent_depth_increased\n"
            )

            for te_name in sorted(rna_data):
                rna_lfc, rna_padj, te_fam, te_cls = rna_data[te_name]
                row_parts = depth_rows.get(te_name)

                if row_parts:
                    # cols: te_name class n_loci wt_norm | samp_norm samp_lfc samp_flag ...
                    cls    = row_parts[1]
                    n_loci = row_parts[2]
                    wt_nd  = row_parts[3]
                    samp_cols = row_parts[4:]   # all per-sample columns

                    flags = [
                        samp_cols[i*3 + 2]
                        for i in range(len(members))
                        if i*3 + 2 < len(samp_cols)
                    ]
                    any_inc        = any(f == "INCREASED" for f in flags)
                    consistent_inc = all(f == "INCREASED" for f in flags) \
                                     and len(flags) > 0
                else:
                    # TE was significant in RNA but had zero depth in WGS
                    # (shouldn't happen often — worth flagging)
                    cls    = te_cls
                    n_loci = "0"
                    wt_nd  = "0"
                    samp_cols      = ["NA\tNA\tNA"] * len(members)
                    any_inc        = False
                    consistent_inc = False

                out.write("\t".join([
                    group, label, te_name, te_fam, te_cls,
                    f"{rna_lfc:.4f}", f"{rna_padj:.4e}",
                    n_loci, wt_nd,
                    "\t".join(samp_cols),
                    str(any_inc), str(consistent_inc),
                ]) + "\n")


rule merge_all_groups:
    """
    Concatenate all 8 group summaries into one master file.
    This is the main deliverable — one table with RNA + DNA evidence
    for every upregulated TE family across all mutants.
    """
    input:
        expand("results/crossref/{group}_summary.tsv", group=GROUPS),
    output:
        "results/crossref/all_groups_summary.tsv",
    threads: 1
    resources:
        mem_mb  = 2000,
        runtime = 5,
    run:
        header_written = False
        with open(output[0], "w") as out:
            for fp in sorted(input):
                with open(fp) as f:
                    for line in f:
                        if line.startswith("#"):
                            continue
                        if line.startswith("group"):
                            if not header_written:
                                out.write(line)
                                header_written = True
                        else:
                            out.write(line)
