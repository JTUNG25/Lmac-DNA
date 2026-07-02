#!/usr/bin/env python3
# =============================================================================
# Logic:
#   1. Read significant_TEs.csv (from R script) → candidate TE family list
#      Only upregulated TEs: log2FC > 1, padj < 0.05
#   2. Filter JN3.te.bed to candidate families only
#   3. fastp  → trim WGS reads
#   4. BWA-MEM → align all samples to plain JN3.fasta
#   5. bedtools coverage → mean read depth per TE locus
#   6. samtools depth → mean depth over actin locus (single-copy reference)
#   7. Collapse loci → per-family actin-normalized depth per sample
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
ACTIN_BED = "data/genome/actin_reference.bed"  # single-copy reference locus

# ── samples ───────────────────────────────────────────────────────────────────
WGS_WT      = "D5"
WGS_MUTANTS = sorted([
    p.name.replace("_R1.fastq.gz", "")
    for p in Path(RAW_DIR).glob("*_R1.fastq.gz")
    if p.name.replace("_R1.fastq.gz", "") != WGS_WT
])
ALL_WGS = WGS_MUTANTS + [WGS_WT]

# ── mutant → R-script label mapping ──────────────────────────────────────────
WGS_TO_LABEL = {
    "A1-1":  "Δago1",  "A1-2":  "Δago1",  "A1-3":  "Δago1",
    "A2-1":  "Δago2",  "A2-2":  "Δago2",  "A2-3":  "Δago2",
    "A3-1":  "Δago3",
    "A13-1": "Δago13", "A13-2": "Δago13",
    "D1":    "Δdcl1",
    "D2-2":  "Δdcl2",  "D2-3":  "Δdcl2",
    "R1-1":  "Δrdrp1",
    "R2-2":  "Δrdrp2", "R2-3":  "Δrdrp2",
    "R2-4":  "Δrdrp2", "R2-5":  "Δrdrp2",
    "R3-1":  "Δrdrp3", "R3-2":  "Δrdrp3", "R3-3":  "Δrdrp3",
    "R12-1": "Δrdrp12","R12-2": "Δrdrp12","R12-3": "Δrdrp12",
}

# ── mutant groups ────────────────────────────────────────────────────────────
GROUPS = {
    "ago1":  {"label": "Δago1",   "samples": ["A1-1", "A1-2", "A1-3"]},
    "ago2":  {"label": "Δago2",   "samples": ["A2-1", "A2-2", "A2-3"]},
    "ago3":  {"label": "Δago3",   "samples": ["A3-1"]},
    "ago13": {"label": "Δago13",  "samples": ["A13-1", "A13-2"]},
    "dcl1":  {"label": "Δdcl1",   "samples": ["D1"]},
    "dcl2":  {"label": "Δdcl2",   "samples": ["D2-2", "D2-3"]},
    "rdrp1": {"label": "Δrdrp1",  "samples": ["R1-1"]},
    "rdrp2": {"label": "Δrdrp2",  "samples": ["R2-2", "R2-3", "R2-4", "R2-5"]},
    "rdrp3": {"label": "Δrdrp3",  "samples": ["R3-1","R3-2", "R3-3"]},
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
        # candidate BED
        "results/te_depth/candidate_te_loci.bed",
        # per-sample depth
        expand("results/te_depth/{sample}.locus_depth.bed",   sample=ALL_WGS),
        expand("results/te_depth/{sample}.actin_depth.txt",   sample=ALL_WGS),
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
        mem_mb  = 8000,
        runtime = 15,
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output}"


rule coverage_summary:
    input:
        expand("results/qc/{sample}.flagstat.txt", sample=ALL_WGS),
    output:
        "results/qc/coverage_summary.txt",
    threads: 1
    resources:
        mem_mb  = 8000,
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
                note   = "WARNING: low read count" if total < 1_000_000 else ""
                out.write(f"{sample}\t{total}\t{mapped}\t{pct}%\t{note}\n")


# =============================================================================
# SECTION 3 — BUILD CANDIDATE TE LOCI BED
# =============================================================================

rule make_candidate_bed:
    input:
        sig_csv = SIG_TE,
        te_bed  = TE_BED,
        fai     = GENOME + ".fai",
    output:
        bed      = "results/te_depth/candidate_te_loci.bed",
        unsorted = temp("results/te_depth/candidate_te_loci.unsorted.bed"),
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 30,
    run:
        candidate_names = set()
        with open(input.sig_csv) as f:
            f.readline()
            for line in f:
                parts = line.strip().split(",")
                if len(parts) < 2:
                    continue
                feature = parts[1].strip('"')
                te_name = feature.split(":")[0]
                candidate_names.add(te_name)

        print(f"  {len(candidate_names)} candidate TE families from R script")

        kept = 0
        with open(input.te_bed) as f_in, open(output.unsorted, "w") as f_out:
            for line in f_in:
                if line.startswith("#") or not line.strip():
                    continue
                col4    = line.strip().split("\t")[3]
                te_name = col4.split(";")[1] if ";" in col4 else col4
                if te_name in candidate_names:
                    f_out.write(line)
                    kept += 1

        print(f"  {kept} TE loci written to unsorted BED")

        shell(
            "module load bedtools/2.31.1-gcc-13.3.0 && "
            "bedtools sort "
            "-i {output.unsorted} "
            "-faidx {input.fai} "
            "> {output.bed}"
        )

rule te_locus_depth:
    input:
        bam = "results/bwa/{sample}.sorted.bam",
        bai = "results/bwa/{sample}.sorted.bam.bai",
        bed = "results/te_depth/candidate_te_loci.bed",
        fai = GENOME + ".fai",
    output:
        "results/te_depth/{sample}.locus_depth.bed",
    container: bedtools
    threads: 2
    resources:
        mem_mb  = 128000,
        runtime = 60,
    shell:
        """
        bedtools coverage \
            -a {input.bed} \
            -b {input.bam} \
            -mean \
            -sorted \
            -g {input.fai} \
        > {output}
        """


# ═══════════════════════════════════════════════════════════════════════════
# NEW RULE: ACTIN DEPTH EXTRACTION
# Extract mean depth over the actin reference locus (Lmb_jn3_12354)
# Single-copy gene used as normalization standard across all samples
# ═══════════════════════════════════════════════════════════════════════════

rule actin_depth:
    """
    Mean read depth over actin locus (Lmb_jn3_12354).
    This is the single-copy reference for normalization.
    bedtools coverage -mean outputs one value per interval.
    """
    input:
        bam = "results/bwa/{sample}.sorted.bam",
        bai = "results/bwa/{sample}.sorted.bam.bai",
        bed = ACTIN_BED,
        fai = GENOME + ".fai",
    output:
        "results/te_depth/{sample}.actin_depth.txt",
    container: bedtools
    threads: 2
    resources:
        mem_mb  = 16000,
        runtime = 15,
    shell:
        """
        bedtools coverage \
            -a {input.bed} \
            -b {input.bam} \
            -mean \
            -sorted \
            -g {input.fai} | \
        awk '{{print $NF}}' > {output}
        """


rule family_depth:
    """
    Collapse per-locus depth → per-family actin-normalized depth.

    For each TE family, average depth across all loci, then divide by
    actin depth (single-copy reference) to normalize for batch sequencing
    depth variation.

    norm_depth = mean_locus_depth / actin_depth
    """
    input:
        depth  = "results/te_depth/{sample}.locus_depth.bed",
        actin  = "results/te_depth/{sample}.actin_depth.txt",
        te_gtf = TE_GTF,
    output:
        "results/te_depth/{sample}.family_depth.tsv",
    threads: 1
    resources:
        mem_mb  = 16000,
        runtime = 10,
    run:
        import re, collections

        with open(input.actin) as f:
            actin_depth = float(f.read().strip() or 0)

        if actin_depth <= 0:
            print(f"WARNING: actin depth is {actin_depth} — check BAM file")
            actin_depth = 1.0  # fallback to avoid division by zero

        # family → class from GTF
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
        fam_depths = collections.defaultdict(list)
        n_loci     = collections.defaultdict(int)
        with open(input.depth) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                col4      = parts[3]
                te_name   = col4.split(";")[1] if ";" in col4 else col4
                locus_dep = float(parts[6])
                fam_depths[te_name].append(locus_dep)
                n_loci[te_name] += 1

        with open(output[0], "w") as out:
            out.write(
                "# Actin-normalised TE family depth\n"
                "# norm_depth = mean_locus_depth / actin_depth\n"
                "te_name\tclass\tsample\tn_loci\t"
                "mean_raw_depth\tactin_depth\tnorm_depth\n"
            )
            for te_name in sorted(fam_depths):
                depths = fam_depths[te_name]
                mean_d = sum(depths) / len(depths)
                norm_d = mean_d / actin_depth if actin_depth > 0 else 0
                cls    = fam_class.get(te_name, "Unknown")
                out.write(
                    f"{te_name}\t{cls}\t{wildcards.sample}\t"
                    f"{n_loci[te_name]}\t{mean_d:.4f}\t{actin_depth:.4f}\t{norm_d:.4f}\n"
                )


# =============================================================================
# SECTION 5 — COMPARE MUTANT DEPTH VS WT
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
        mem_mb  = 8000,
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
                "# ACTIN-NORMALIZED depth comparison\n"
                "# log2FC_depth = log2(mutant_norm / wt_norm)\n"
                "# All depths normalized by actin (Lmb_jn3_12354)\n"
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
# =============================================================================

rule crossref_group:
    input:
        sig_csv   = SIG_TE,
        depth_tsv = "results/te_depth/{group}_vs_wt.tsv",
    output:
        "results/crossref/{group}_summary.tsv",
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 5,
    run:
        import math

        group   = wildcards.group
        label   = GROUPS[group]["label"]
        members = GROUPS[group]["samples"]

        # RNA data for this mutant label
        rna_data = {}
        with open(input.sig_csv) as f:
            f.readline()
            for line in f:
                parts   = line.strip().split(",")
                if len(parts) < 6:
                    continue
                mutant  = parts[0].strip('"')
                feature = parts[1].strip('"')
                te_name = parts[2].strip('"')
                te_fam  = parts[3].strip('"')
                te_cls  = parts[4].strip('"')
                try:
                    lfc  = float(parts[5])
                    padj = float(parts[8])
                except ValueError:
                    continue
                rna_data[te_name] = (lfc, padj, te_fam, te_cls)

        # Depth results
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
                    depth_rows[parts[0]] = parts

        # Write output
        samp_hdr = "\t".join(
            [f"{s}_norm_depth\t{s}_log2FC_depth\t{s}_depth_flag"
             for s in members]
        )
        with open(output[0], "w") as out:
            out.write(
                "# Cross-reference: RNA upregulation (TEtranscripts) vs "
                "WGS read depth (ACTIN-NORMALIZED)\n"
                "# All depths normalized by actin (Lmb_jn3_12354)\n"
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
                    cls    = row_parts[1]
                    n_loci = row_parts[2]
                    wt_nd  = row_parts[3]
                    samp_cols = row_parts[4:]

                    flags = [
                        samp_cols[i*3 + 2]
                        for i in range(len(members))
                        if i*3 + 2 < len(samp_cols)
                    ]
                    any_inc        = any(f == "INCREASED" for f in flags)
                    consistent_inc = all(f == "INCREASED" for f in flags) \
                                     and len(flags) > 0
                else:
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
    input:
        expand("results/crossref/{group}_summary.tsv", group=GROUPS),
    output:
        "results/crossref/all_groups_summary.tsv",
    threads: 1
    resources:
        mem_mb  = 8000,
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
