#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
BAMDIR = "results/3_mismatch/2_marked_new_reads"
OUTDIR = "results/4_expression"

rule all:
    input:
        expand(OUTDIR + "/1_fpkm/{run_cell}.tsv", run_cell=RUN_CELLS),

rule calculate_fpkm:
    input:
        bam = BAMDIR + "/{run}/{cell}.bam",
        bed = config["TRANSCRIPT_BED_GZ"],
        tsv = config["TRANSCRIPT_INFO_TSV"]
    output:
        txt = OUTDIR + "/1_fpkm/{run}/{cell}.tsv"
    log:
        OUTDIR + "/1_fpkm/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        calculate_fpkm.py --threads {threads} --annotation {input.tsv} \
            --stranded TAG --strand-tag ST --new \
            {input.bam} {input.bed} {output.txt} &> {log}
        """
