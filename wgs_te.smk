#!/usr/bin/env python3
# =============================================================================
# wgs_te.smk  —  WGS TE copy-number & insertion screen
# Leptosphaeria maculans JN3
#
# Pipeline:
#   0. fastp          — trim raw WGS reads
#   1. BWA-MEM        — align to JN3 reference (single consistent reference)
#   2. samtools       — sort, index, QC flagstat
#   3. TE depth       — per-family normalised depth vs WT (exploratory, no stats)
#   4. McClintock 2   — non-reference insertion detection
#   5. Cross-ref      — intersect with TEtranscripts sig families
#
# NOTE: Single WT sample, single mutant samples — no replicates.
#       Depth comparisons are EXPLORATORY / SCREENING only.
#       Flagged candidates should be validated with biological replicates.
# =============================================================================

import os
import re
from pathlib import Path

# ── containers ────────────────────────────────────────────────────────────────
fastp      = "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_3"
bwa        = "docker://quay.io/biocontainers/bwa:0.7.18--he4a0461_0"
samtools   = "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
bedtools   = "docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"
mcclintock = "docker://docker.io/bergmanlab/mcclintock:2.0.3"

# ── paths ─────────────────────────────────────────────────────────────────────
RAW_DIR      = "/QRISdata/Q9141/lmac_dna/DNAseq"
GENOME       = "data/genome/JN3.fasta"
TE_GTF       = "data/genome/JN3.te.gtf"
TE_BED       = "data/genome/JN3.te.bed"
TE_CONSENSUS = "data/repeatmodeler/JN3-families.fa"

# ── sample discovery from raw FASTQs ─────────────────────────────────────────
# Pattern: {sample}_{flowcell}_{barcodes}_{lane}_R1.fastq.gz
# We extract the sample name (everything before the first _23)
_raw_r1 = list(Path(RAW_DIR).glob("*_R1.fastq.gz"))

def _extract_sample(path):
    """Return sample name from raw FASTQ filename."""
    return re.split(r"_\d{2}[A-Z]", path.name)[0]

# build sample -> {R1, R2} mapping
SAMPLE_FILES = {}
for r1 in _raw_r1:
    sname = _extract_sample(r1)
    r2 = Path(str(r1).replace("_R1.fastq.gz", "_R2.fastq.gz"))
    if r2.exists():
        SAMPLE_FILES[sname] = {"R1": str(r1), "R2": str(r2)}

# ── sample lists ──────────────────────────────────────────────────────────────
# WT is the D5 sample (single replicate — screening only)
WGS_WT = "D5"

WGS_MUTANTS = sorted([
    s for s in SAMPLE_FILES
])

ALL_WGS = WGS_MUTANTS + [WGS_WT]

# mutant groups — mirrors TEtranscripts batches
MUTANT_GROUPS = {
    "A1":  [s for s in WGS_MUTANTS if s.startswith("A1-")],
    "A3":  [s for s in WGS_MUTANTS if s.startswith("A3-")],
    "A13": [s for s in WGS_MUTANTS if s.startswith("A13-")],
    "D2":  [s for s in WGS_MUTANTS if s.startswith("D2-")],
    "R1":  [s for s in WGS_MUTANTS if s.startswith("R1-")],
    "R2":  [s for s in WGS_MUTANTS if s.startswith("R2-")],
    "R12": [s for s in WGS_MUTANTS if s.startswith("R12-")],
}
# drop empty groups (e.g. if some samples missing from DNAseq dir)
MUTANT_GROUPS = {k: v for k, v in MUTANT_GROUPS.items() if v}

# ── TEtranscripts results (from RNA project) ──────────────────────────────────
RNA_BASE = "/QRISdata/Q9141/lmac_rna/results/tetranscripts"
TETRANSCRIPTS_RESULTS = {
    "A1":  f"{RNA_BASE}/batch2a/lepto_batch2a_DESeq_TE_results.txt",
    "A3":  f"{RNA_BASE}/batch2b/lepto_batch2b_DESeq_TE_results.txt",
    "A13": f"{RNA_BASE}/batch4a/lepto_batch4a_DESeq_TE_results.txt",
    "D2":  f"{RNA_BASE}/batch1/lepto_batch1_DESeq_TE_results.txt",
    "R1":  f"{RNA_BASE}/batch3a/lepto_batch3a_DESeq_TE_results.txt",
    "R2":  f"{RNA_BASE}/batch3b/lepto_batch3b_DESeq_TE_results.txt",
    "R12": f"{RNA_BASE}/batch4b/lepto_batch4b_DESeq_TE_results.txt",
}
TETRANSCRIPTS_RESULTS = {
    k: v for k, v in TETRANSCRIPTS_RESULTS.items() if k in MUTANT_GROUPS
}

# ── helpers ───────────────────────────────────────────────────────────────────
def raw_r1(wc): return SAMPLE_FILES[wc.sample]["R1"]
def raw_r2(wc): return SAMPLE_FILES[wc.sample]["R2"]

rule all:
    input:
        # Section 0+1+2 — trim, align, QC
        expand("results/fastp/{sample}_R1.fastq.gz",   sample=ALL_WGS),
        expand("results/bwa/{sample}.sorted.bam",      sample=ALL_WGS),
        expand("results/bwa/{sample}.sorted.bam.bai",  sample=ALL_WGS),
        expand("results/qc/{sample}.flagstat.txt",     sample=ALL_WGS),
        "results/qc/coverage_summary.txt",
        # Section 3 — TE depth
        expand("results/te_depth/{sample}.family_depth.tsv", sample=ALL_WGS),
        expand("results/te_depth/{group}_depth_vs_wt.tsv",   group=MUTANT_GROUPS),
        # Section 4 — McClintock
        expand("results/mcclintock/{sample}/results/summary/", sample=ALL_WGS),
        # Section 5 — cross-reference
        expand("results/crossref/{group}_candidates.tsv", group=MUTANT_GROUPS),
        "results/crossref/all_groups_summary.tsv",

rule fastp_trim:
    input:
        r1 = raw_r1,
        r2 = raw_r2,
    output:
        r1      = "results/fastp/{sample}_R1.fastq.gz",
        r2      = "results/fastp/{sample}_R2.fastq.gz",
        html    = "results/fastp/{sample}_fastp.html",
        json    = "results/fastp/{sample}_fastp.json",
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
# All samples aligned to the same plain JN3.fasta — no T-DNA reference.
# This is essential so depth comparisons between mutants and WT are valid.
# =============================================================================

rule bwa_index:
    """Build BWA index — only needed once, fast on 45 Mb genome."""
    input:  GENOME,
    output:
        multiext(GENOME, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    container: bwa
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 30,
    shell:
        "bwa index {input}"


rule bwa_mem_align:
    """
    BWA-MEM flags chosen for TE analysis:
      -a  report all alignments for multi-mappers
          (McClintock needs split-read evidence at TE junctions)
      -M  mark shorter split hits as secondary
          (compatibility with downstream tools)
    Read group added so BAMs are traceable.
    """
    input:
        r1    = "results/fastp/{sample}_R1.fastq.gz",
        r2    = "results/fastp/{sample}_R2.fastq.gz",
        index = multiext(GENOME, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam = "results/bwa/{sample}.sorted.bam",
        bai = "results/bwa/{sample}.sorted.bam.bai",
    container: samtools   # samtools sif used for sort/index; bwa called via bwa sif below
    threads: 16
    resources:
        mem_mb  = 64000,
        runtime = 240,
    params:
        genome  = GENOME,
        rg      = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA",
        tmp     = "/scratch/temp/$SLURM_JOB_ID/bwa_{sample}",
        bwa_sif = bwa,
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
        samtools sort \
                -@ {threads} \
                -T {params.tmp}/sort \
                -o {output.bam}

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
    output: "results/qc/{sample}.flagstat.txt",
    container: samtools
    threads: 2
    resources:
        mem_mb  = 4000,
        runtime = 15,
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output}"


rule coverage_summary:
    """
    Aggregate flagstat across all samples.
    Low-coverage samples (< 1 M reads) are flagged with a WARNING.
    This is where A1-3 (and any other problematic sample) will show up.
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
                fpath = f"results/qc/{sample}.flagstat.txt"
                with open(fpath) as f:
                    lines = f.readlines()
                total  = int(lines[0].split()[0])
                mapped = int([l for l in lines if "mapped (" in l][0].split()[0])
                pct    = f"{100*mapped/total:.1f}" if total > 0 else "0"
                note   = "WARNING: very low read count — exclude from analysis" \
                         if total < 1_000_000 else ""
                out.write(f"{sample}\t{total}\t{mapped}\t{pct}%\t{note}\n")


# =============================================================================
# SECTION 3 — TE READ DEPTH (family-level copy-number proxy)
#
# Strategy:
#   a) bedtools coverage -mean → mean depth per TE locus in JN3.te.bed
#   b) samtools coverage       → genome-wide mean depth (for normalisation)
#   c) collapse loci → family mean normalised depth
#   d) compare each mutant to WT → log2 fold-change
#
# IMPORTANT: No replicates → no statistical test.
# FC >= 1.5 (log2FC >= 0.585) used as a screening threshold only.
# Results labelled EXPLORATORY in output headers.
# =============================================================================

rule te_locus_depth:
    input:
        bam = "results/bwa/{sample}.sorted.bam",
        bai = "results/bwa/{sample}.sorted.bam.bai",
        bed = TE_BED,
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
            -sorted \
        > {output}
        """


rule genome_mean_depth:
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
    Collapse per-locus depth to per-family normalised depth.
    BED col4 format (from make_te_gtf in mrna_te.smk): locus_id;family_name
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

        # family -> class map from GTF
        fam_class = {}
        with open(input.te_gtf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                mg = re.search(r'gene_id "([^"]+)"', line)
                mc = re.search(r'class_id "([^"]+)"', line)
                if mg and mc:
                    fam_class[mg.group(1)] = mc.group(1)

        fam_depths = collections.defaultdict(list)
        with open(input.depth) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                locus_fam = parts[3]
                raw_depth = float(parts[6])
                fam = locus_fam.split(";")[1] if ";" in locus_fam else locus_fam
                fam_depths[fam].append(raw_depth)

        with open(output[0], "w") as out:
            out.write(
                "# EXPLORATORY SCREENING — single WT, no replicates, no statistical test\n"
                "family\tclass\tsample\tmean_raw_depth\tgenome_mean_depth\tnorm_depth\n"
            )
            for fam in sorted(fam_depths):
                depths  = fam_depths[fam]
                mean_d  = sum(depths) / len(depths)
                norm_d  = mean_d / gd if gd > 0 else 0
                cls     = fam_class.get(fam, "Unknown")
                out.write(
                    f"{fam}\t{cls}\t{wildcards.sample}\t"
                    f"{mean_d:.4f}\t{gd:.4f}\t{norm_d:.4f}\n"
                )


rule depth_vs_wt:
    """
    Per-group: compare each mutant's normalised depth to WT.
    log2FC = log2(mutant_norm / wt_norm)
    Threshold >= log2(1.5) flagged as INCREASED — screening only.
    """
    input:
        mutant_depths = lambda wc: expand(
            "results/te_depth/{sample}.family_depth.tsv",
            sample=MUTANT_GROUPS[wc.group]
        ),
        wt_depth = f"results/te_depth/{WGS_WT}.family_depth.tsv",
    output:
        "results/te_depth/{group}_depth_vs_wt.tsv",
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 5,
    run:
        import math, collections, re

        group   = wildcards.group
        members = MUTANT_GROUPS[group]

        def parse_depth_file(path):
            """Return {family: (class, norm_depth)} skipping comment lines."""
            result = {}
            with open(path) as f:
                for line in f:
                    if line.startswith("#") or line.startswith("family"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 6:
                        continue
                    fam, cls, samp, _, _, norm = parts
                    result[fam] = (cls, float(norm))
            return result

        wt_data = parse_depth_file(input.wt_depth)

        # each mutant gets its own column in the output
        mutant_data = {}
        for fp in input.mutant_depths:
            sname = re.search(r"te_depth/(.+)\.family_depth", fp).group(1)
            mutant_data[sname] = parse_depth_file(fp)

        all_families = set(wt_data)
        for d in mutant_data.values():
            all_families |= set(d)

        with open(output[0], "w") as out:
            # header
            mut_cols = "\t".join(
                [f"{s}_norm_depth\t{s}_log2FC\t{s}_flag" for s in members]
            )
            out.write(
                "# EXPLORATORY SCREENING — no replicates, no statistical testing\n"
                "# Flag threshold: log2FC >= 0.585 (1.5-fold increase vs single WT)\n"
                f"family\tclass\twt_norm_depth\t{mut_cols}\n"
            )
            for fam in sorted(all_families):
                wt_cls, wt_norm = wt_data.get(fam, ("Unknown", 0.0))
                row = [fam, wt_cls, f"{wt_norm:.4f}"]
                for samp in members:
                    _, mut_norm = mutant_data[samp].get(fam, ("Unknown", 0.0))
                    if wt_norm > 0 and mut_norm > 0:
                        lfc = math.log2(mut_norm / wt_norm)
                    elif mut_norm > 0:
                        lfc = float("inf")
                    else:
                        lfc = float("-inf")
                    flag = "INCREASED" if lfc != float("inf") and lfc >= math.log2(1.5) else \
                           "NEW"       if lfc == float("inf") else ""
                    row += [f"{mut_norm:.4f}", f"{lfc:.4f}", flag]
                out.write("\t".join(row) + "\n")


# =============================================================================
# SECTION 4 — McCLINTOCK 2
# Detects non-reference TE insertions using split-read + read-pair evidence.
# Uses four complementary methods internally:
#   ngs_te_mapper2, RelocaTE2, TEMP2, TEbreak
# BWA-MEM alignments are required here — Bowtie2 alignments would miss
# split reads at TE-genome junctions.
# =============================================================================

rule mcclintock:
    input:
        bam    = "results/bwa/{sample}.sorted.bam",
        bai    = "results/bwa/{sample}.sorted.bam.bai",
        genome = GENOME,
        te_lib = TE_CONSENSUS,
    output:
        summary = directory("results/mcclintock/{sample}/results/summary/"),
    container: mcclintock
    threads: 8
    resources:
        mem_mb  = 32000,
        runtime = 480,
    params:
        outdir = "results/mcclintock/{sample}",
        sample = "{sample}",
    shell:
        """
        mkdir -p {params.outdir}

        mcclintock.py \
            --bam         {input.bam} \
            --reference   {input.genome} \
            --consensus   {input.te_lib} \
            --sample_name {params.sample} \
            --out         {params.outdir} \
            --proc        {threads} \
            --method      ngs_te_mapper2 relocate2 temp2 tebreak \
            --clean
        """


# =============================================================================
# SECTION 5 — CROSS-REFERENCE
#
# For each mutant group, intersect:
#   - sig upregulated TE families (TEtranscripts, padj<0.05, log2FC>0)
#   - depth fold-change vs WT (Section 3)
#   - McClintock insertion calls (Section 4)
#
# Name matching:
#   TEtranscripts: gene_id = RepeatModeler name e.g. "Copia-1_LM"
#   McClintock:    consensus header e.g. "Copia-1_LM#LTR/Copia"
#   → strip after "#", then exact match; fallback to 8-char prefix match
#   Both original names kept in output so mismatches are auditable.
#
# Evidence levels:
#   RNA_only              → transcriptional derepression, no DNA change
#   RNA+depth             → copy number elevated, no mapped insertion site
#   RNA+insertion         → new insertion site(s) found
#   RNA+depth+insertion   → strongest evidence of active transposition
# =============================================================================

rule crossref:
    input:
        tet      = lambda wc: TETRANSCRIPTS_RESULTS[wc.group],
        depth    = "results/te_depth/{group}_depth_vs_wt.tsv",
        mc_dirs  = lambda wc: expand(
            "results/mcclintock/{sample}/results/summary/",
            sample=MUTANT_GROUPS[wc.group]
        ),
    output:
        "results/crossref/{group}_candidates.tsv",
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 10,
    run:
        import glob, re, math, collections

        group   = wildcards.group
        members = MUTANT_GROUPS[group]

        # ── 1. Sig upregulated families from TEtranscripts ────────────────
        # Format: gene_id  baseMean  log2FC  lfcSE  stat  pvalue  padj
        sig = {}   # name -> (log2FC, padj)
        with open(input.tet) as f:
            f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                name = parts[0].strip('"')
                try:
                    lfc  = float(parts[2])
                    padj = float(parts[6])
                except ValueError:
                    continue
                if padj < 0.05 and lfc > 0:
                    sig[name] = (lfc, padj)

        # ── 2. Depth table: family -> {sample -> (norm_depth, log2FC, flag)} ─
        depth_by_fam = collections.defaultdict(dict)
        depth_class  = {}
        wt_depth_by_fam = {}
        with open(input.depth) as f:
            header = None
            for line in f:
                if line.startswith("#"):
                    continue
                if line.startswith("family"):
                    header = line.strip().split("\t")
                    continue
                parts = line.strip().split("\t")
                fam   = parts[0]
                cls   = parts[1]
                depth_class[fam] = cls
                wt_depth_by_fam[fam] = float(parts[2])
                # columns: family class wt_norm | sample_norm sample_lfc sample_flag ...
                for i, samp in enumerate(members):
                    col = 3 + i * 3
                    if col + 2 < len(parts):
                        depth_by_fam[fam][samp] = {
                            "norm":  float(parts[col]),
                            "log2FC": parts[col+1],
                            "flag":   parts[col+2],
                        }

        # ── 3. McClintock insertions: te_name -> {sample -> count} ────────
        mc_counts = collections.defaultdict(lambda: collections.defaultdict(int))
        for samp in members:
            pattern = f"results/mcclintock/{samp}/results/summary/*.bed"
            for bed in glob.glob(pattern):
                with open(bed) as f:
                    for line in f:
                        if line.startswith("#") or not line.strip():
                            continue
                        p = line.strip().split("\t")
                        if len(p) >= 4:
                            mc_counts[p[3]][samp] += 1

        # ── 4. Name matching ──────────────────────────────────────────────
        def norm_name(n):
            return re.sub(r"#.*", "", n).strip().lower()

        mc_norm_map = {norm_name(k): k for k in mc_counts}

        def find_mc(tet_name):
            n = norm_name(tet_name)
            if n in mc_norm_map:
                return mc_norm_map[n]
            for mc_n, mc_orig in mc_norm_map.items():
                if len(n) >= 8 and len(mc_n) >= 8:
                    if mc_n.startswith(n[:8]) or n.startswith(mc_n[:8]):
                        return mc_orig
            return None

        # ── 5. Write output ───────────────────────────────────────────────
        with open(output[0], "w") as out:
            out.write(
                "# EXPLORATORY SCREENING — single WT replicate, no statistical testing on depth\n"
                "# Candidates should be validated with biological replicates\n"
            )
            # per-sample depth/insertion columns
            samp_depth_hdr = "\t".join(
                [f"{s}_log2FC_depth\t{s}_depth_flag\t{s}_mc_insertions"
                 for s in members]
            )
            out.write(
                f"group\ttet_family\tmc_family\tclass\t"
                f"tet_log2FC\ttet_padj\twt_norm_depth\t"
                f"{samp_depth_hdr}\t"
                f"any_depth_increased\ttotal_mc_insertions\t"
                f"mc_samples_with_hits\tevidence_level\n"
            )

            for fam, (tet_lfc, tet_padj) in sorted(sig.items()):
                cls     = depth_class.get(fam, "Unknown")
                wt_nd   = wt_depth_by_fam.get(fam, 0.0)
                mc_key  = find_mc(fam)

                any_increased  = False
                total_mc_ins   = 0
                mc_samp_hits   = 0
                samp_cols      = []

                for samp in members:
                    dep   = depth_by_fam[fam].get(samp, {})
                    dlfc  = dep.get("log2FC", "NA")
                    dflag = dep.get("flag",   "")
                    if dflag == "INCREASED":
                        any_increased = True

                    ins = 0
                    if mc_key:
                        ins = mc_counts[mc_key].get(samp, 0)
                        total_mc_ins += ins
                        if ins > 0:
                            mc_samp_hits += 1

                    samp_cols.append(f"{dlfc}\t{dflag}\t{ins}")

                has_depth = any_increased
                has_ins   = total_mc_ins > 0
                if has_depth and has_ins:
                    evidence = "RNA+depth+insertion"
                elif has_depth:
                    evidence = "RNA+depth"
                elif has_ins:
                    evidence = "RNA+insertion"
                else:
                    evidence = "RNA_only"

                out.write("\t".join([
                    group, fam, mc_key or "NA", cls,
                    f"{tet_lfc:.4f}", f"{tet_padj:.4e}",
                    f"{wt_nd:.4f}",
                    "\t".join(samp_cols),
                    str(any_increased), str(total_mc_ins),
                    str(mc_samp_hits), evidence,
                ]) + "\n")


rule merge_crossref:
    """Master summary across all groups — one file for your supervisors."""
    input:
        expand("results/crossref/{group}_candidates.tsv", group=MUTANT_GROUPS),
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
