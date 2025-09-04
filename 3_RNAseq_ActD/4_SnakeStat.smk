#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
BAMDIR = "results/02_mapping/04_rmdup"
OUTDIR = "results/04_stat"

rule all:
    input:
        expand(OUTDIR + "/01_read_location/{sample}_summary.tsv", sample=SAMPLES),

rule stat_read_location:
    input:
        bam = BAMDIR + "/{sample}.human.bam",
        bed = config["HUMAN_TRANSCRIPT_BED"]
    output:
        bed1 = temp(OUTDIR + "/01_read_location/{sample}.bed"),
        bed2 = OUTDIR + "/01_read_location/{sample}.bed.gz",
        tsv = OUTDIR + "/01_read_location/{sample}_summary.tsv"
    log:
        OUTDIR + "/01_read_location/{sample}.log"
    threads:
        8
    shell:
        """
        stat_read_location.py -b {input.bam} -g {input.bed} -t {threads} -o {output.bed1} -s {output.tsv} &> {log}
        bgzip -c {output.bed1} > {output.bed2}
        """