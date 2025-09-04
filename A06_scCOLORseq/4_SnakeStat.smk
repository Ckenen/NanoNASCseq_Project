#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
BAMDIR = "results/03_mapping/03_bam_filtered"
OUTDIR = "results/04_stat"

rule all:
    input:
        expand(OUTDIR + "/01_read_location/{sample}_summary.tsv", sample=SAMPLES)

rule stat_read_location:
    input:
        bam = BAMDIR + "/{sample}.bam",
        bed = config["TRANSCRIPT_BED_GZ"]
    output:
        bed1 = temp(OUTDIR + "/01_read_location/{sample}.bed"),
        bed2 = OUTDIR + "/01_read_location/{sample}.bed.gz",
        tsv = OUTDIR + "/01_read_location/{sample}_summary.tsv",
    log:
        OUTDIR + "/01_read_location/{sample}.log"
    shell:
        """
        stat_read_location.py -l {input.bam} -g {input.bed} -o {output.bed1} -s {output.tsv} &> {log}
        bgzip -c {output.bed1} > {output.bed2}
        """