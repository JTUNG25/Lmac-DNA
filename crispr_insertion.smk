#!/usr/bin/env python3
# =============================================================================
# Workflow:
#   - Bowtie2 align (using pre-trimmed reads) → samtools sort/index
#   - Extract discordant, soft-clipped, unmapped reads
#   - Generate insertion hotspot map
#   - Produce QC report and candidate insertion sites table
# =============================================================================

import os
import re
from pathlib import Path
from collections import defaultdict

# ── paths ─────────────────────────────────────────────────────────────────────
FASTP_DIR = "results/fastp"
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
        expand("{genome}.1.bt2", genome=ALL_REFERENCES),
        expand("results/bowtie2/{sample}.sorted.bam", sample=ALL_SAMPLES),
        expand("results/bowtie2/{sample}.sorted.bam.bai", sample=ALL_SAMPLES),
        expand("results/qc/{sample}.flagstat.txt", sample=ALL_SAMPLES),
        "results/qc/mapping_summary.txt",
        expand("results/reads/{sample}.discordant.bam", sample=ALL_SAMPLES),
        expand("results/reads/{sample}.softclipped.bam", sample=ALL_SAMPLES),
        expand("results/reads/{sample}.unmapped_mate.bam", sample=ALL_SAMPLES),
        expand("results/hotspots/{sample}.insertion_hotspots.bed", sample=ALL_SAMPLES),
        expand("results/summary/{sample}.insertion_candidates.tsv", sample=ALL_SAMPLES),
        "results/summary/all_samples_candidates.tsv",


# =============================================================================
# SECTION 1 — BOWTIE2 ALIGNMENT
# =============================================================================


rule bowtie2_index:
    input:
        "{genome}",
    output:
        multiext(
            "{genome}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    threads: 8
    resources:
        mem_mb=32000,
        runtime=60,
    params:
        index_name="{genome}",
    shell:
        """
        module load bowtie2/2.5.1-gcc-12.3.0
        bowtie2-build --threads {threads} {input} {params.index_name}
        """


rule bowtie2_align:
    input:
        r1=r1,
        r2=r2,
        ref=get_ref,
        idx1=lambda wc: get_ref(wc) + ".1.bt2",
        idx2=lambda wc: get_ref(wc) + ".2.bt2",
    output:
        bam="results/bowtie2/{sample}.sorted.bam",
        bai="results/bowtie2/{sample}.sorted.bam.bai",
    threads: 16
    resources:
        mem_mb=64000,
        runtime=240,
    shell:
        """
        module load bowtie2/2.5.1-gcc-12.3.0
        module load samtools/1.18-gcc-12.3.0

        mkdir -p tmp/bowtie2_{wildcards.sample}
        TMPDIR=tmp/bowtie2_{wildcards.sample}

        bowtie2 -x {input.ref} \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} \
            -k 5 -X 1000 --very-sensitive \
            -S $TMPDIR/{wildcards.sample}.sam

        samtools sort -@ {threads} -T $TMPDIR/sort \
            $TMPDIR/{wildcards.sample}.sam \
            -o {output.bam}

        samtools index -@ {threads} {output.bam}

        rm -rf $TMPDIR
        """


# =============================================================================
# SECTION 2 — QC AND MAPPING STATISTICS
# =============================================================================


rule flagstat:
    input:
        bam="results/bowtie2/{sample}.sorted.bam",
        bai="results/bowtie2/{sample}.sorted.bam.bai",
    output:
        "results/qc/{sample}.flagstat.txt",
    threads: 2
    resources:
        mem_mb=8000,
        runtime=15,
    shell:
        """
        module load samtools/1.18-gcc-12.3.0
        samtools flagstat -@ {threads} {input.bam} >{output}
        """


rule mapping_summary:
    input:
        expand("results/qc/{sample}.flagstat.txt", sample=ALL_SAMPLES),
    output:
        "results/qc/mapping_summary.txt",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=5,
    run:
        with open(output[0], "w") as out:
            out.write(
                "sample\ttotal_reads\tmapped_reads\tmapped_pct\t"
                "properly_paired\tproperly_paired_pct\tnote\n"
            )
            for sample in sorted(ALL_SAMPLES):
                with open(f"results/qc/{sample}.flagstat.txt") as f:
                    lines = f.readlines()
                total = int(lines[0].split()[0])
                mapped = int([l for l in lines if "mapped (" in l][0].split()[0])
                properly = int(
                    [l for l in lines if "properly paired" in l][0].split()[0]
                )
                mapped_pct = f"{100*mapped/total:.1f}" if total > 0 else "0"
                properly_pct = f"{100*properly/total:.1f}" if total > 0 else "0"
                note = ""
                if total < 1_000_000:
                    note = "WARNING: low read count"
                elif mapped / total < 0.7:
                    note = "WARNING: low mapping rate"
                out.write(
                    f"{sample}\t{total}\t{mapped}\t{mapped_pct}%\t"
                    f"{properly}\t{properly_pct}%\t{note}\n"
                )


# =============================================================================
# SECTION 3 — EXTRACT READS WITH INSERTION EVIDENCE
# =============================================================================


rule extract_discordant:
    input:
        bam="results/bowtie2/{sample}.sorted.bam",
        bai="results/bowtie2/{sample}.sorted.bam.bai",
    output:
        "results/reads/{sample}.discordant.bam",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=30,
    shell:
        """
        module load samtools/1.18-gcc-12.3.0
        samtools view -b -@ {threads} -F 2 {input.bam} >{output}
        """


rule extract_softclipped:
    input:
        bam="results/bowtie2/{sample}.sorted.bam",
        bai="results/bowtie2/{sample}.sorted.bam.bai",
    output:
        "results/reads/{sample}.softclipped.bam",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=30,
    run:
        import pysam

        inbam = pysam.AlignmentFile(input.bam, "rb")
        outbam = pysam.AlignmentFile(output[0], "wb", template=inbam)
        count = 0
        for read in inbam:
            if read.cigarstring and "S" in read.cigarstring:
                outbam.write(read)
                count += 1
        inbam.close()
        outbam.close()


rule extract_unmapped_mate:
    input:
        bam="results/bowtie2/{sample}.sorted.bam",
        bai="results/bowtie2/{sample}.sorted.bam.bai",
    output:
        "results/reads/{sample}.unmapped_mate.bam",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=30,
    shell:
        """
        module load samtools/1.18-gcc-12.3.0
        samtools view -b -@ {threads} -f 4 -F 8 {input.bam} >{output}
        """


# =============================================================================
# SECTION 4 — IDENTIFY INSERTION HOTSPOTS
# =============================================================================


rule find_hotspots:
    input:
        discordant="results/reads/{sample}.discordant.bam",
        softclipped="results/reads/{sample}.softclipped.bam",
        unmapped="results/reads/{sample}.unmapped_mate.bam",
        fai=lambda wc: get_sample_reference(wc) + ".fai",
    output:
        hotspots="results/hotspots/{sample}.insertion_hotspots.bed",
    threads: 2
    resources:
        mem_mb=16000,
        runtime=30,
    run:
        import pysam
        from collections import defaultdict

        positions = defaultdict(
            lambda: {"discordant": 0, "softclipped": 0, "unmapped": 0}
        )
        try:
            bam = pysam.AlignmentFile(input.discordant, "rb")
            for read in bam:
                if not read.is_unmapped:
                    key = f"{read.reference_name}:{read.reference_start}"
                    positions[key]["discordant"] += 1
            bam.close()
        except:
            pass
        try:
            bam = pysam.AlignmentFile(input.softclipped, "rb")
            for read in bam:
                if not read.is_unmapped:
                    key = f"{read.reference_name}:{read.reference_start}"
                    positions[key]["softclipped"] += 1
            bam.close()
        except:
            pass
        try:
            bam = pysam.AlignmentFile(input.unmapped, "rb")
            for read in bam:
                if (
                    read.next_reference_name
                    and not read.next_reference_name.startswith("*")
                ):
                    key = f"{read.next_reference_name}:{read.next_reference_start}"
                    positions[key]["unmapped"] += 1
            bam.close()
        except:
            pass
        with open(output.hotspots, "w") as out:
            out.write("# Insertion evidence hotspots\n")
            out.write(
                "# chrom\tstart\tend\tevidence_score\tdiscordant_count\t"
                "softclipped_count\tunmapped_mate_count\n"
            )
            for key in sorted(positions):
                chrom, pos = key.split(":")
                pos = int(pos)
                evidence = positions[key]
                score = (
                    evidence["discordant"] * 2
                    + evidence["softclipped"] * 1
                    + evidence["unmapped"] * 1
                )
                if score >= 2:
                    out.write(
                        f"{chrom}\t{pos}\t{pos+1}\t{score}\t"
                        f"{evidence['discordant']}\t"
                        f"{evidence['softclipped']}\t"
                        f"{evidence['unmapped']}\n"
                    )


# =============================================================================
# SECTION 5 — GENERATE INSERTION CANDIDATE SUMMARY
# =============================================================================


rule candidate_summary:
    input:
        hotspots="results/hotspots/{sample}.insertion_hotspots.bed",
        bam="results/bowtie2/{sample}.sorted.bam",
        bai="results/bowtie2/{sample}.sorted.bam.bai",
    output:
        "results/summary/{sample}.insertion_candidates.tsv",
    threads: 2
    resources:
        mem_mb=8000,
        runtime=15,
    run:
        import pysam

        hotspots_list = []
        with open(input.hotspots) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 7:
                    chrom, start, end, score, disc, soft, unmapped = parts[:7]
                    hotspots_list.append(
                        {
                            "chrom": chrom,
                            "pos": int(start),
                            "score": int(score),
                            "discordant": int(disc),
                            "softclipped": int(soft),
                            "unmapped": int(unmapped),
                        }
                    )
        hotspots_list.sort(key=lambda x: x["score"], reverse=True)
        bam = pysam.AlignmentFile(input.bam, "rb")
        with open(output[0], "w") as out:
            out.write(
                "rank\tchrom\tposition\ttotal_evidence_score\t"
                "discordant_reads\tsoftclipped_reads\tunmapped_mate_reads\t"
                "local_coverage\tpriority\n"
            )
            for rank, hs in enumerate(hotspots_list, 1):
                try:
                    coverage = bam.count(hs["chrom"], hs["pos"], hs["pos"] + 1)
                except:
                    coverage = "NA"
                if hs["score"] >= 10 and hs["discordant"] >= 3:
                    priority = "HIGH"
                elif hs["score"] >= 5 and hs["softclipped"] >= 2:
                    priority = "MEDIUM"
                else:
                    priority = "LOW"
                out.write(
                    f"{rank}\t{hs['chrom']}\t{hs['pos']}\t{hs['score']}\t"
                    f"{hs['discordant']}\t{hs['softclipped']}\t"
                    f"{hs['unmapped']}\t{coverage}\t{priority}\n"
                )
        bam.close()


rule merge_all_candidates:
    input:
        expand("results/summary/{sample}.insertion_candidates.tsv", sample=ALL_SAMPLES),
    output:
        "results/summary/all_samples_candidates.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=5,
    run:
        header_written = False
        with open(output[0], "w") as out:
            out.write("sample\t")
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
