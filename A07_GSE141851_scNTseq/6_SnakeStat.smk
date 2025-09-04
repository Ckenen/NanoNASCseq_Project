#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = list(filter(lambda s: "K562" in s, config["SAMPLES"]))
BAMDIR = "results/02_dropseq/12_DetectedErrors"
OUTDIR = "results/06_stat"

rule all:
    input:
        expand(OUTDIR + "/01_sorted_bam/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/02_read_location/{sample}_summary.tsv", sample=SAMPLES),

rule get_bam:
    input:
        bam = BAMDIR + "/{sample}.bam"
    output:
        bam = OUTDIR + "/01_sorted_bam/{sample}.bam"
    threads:
        4
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule stat_read_location:
    input:
        bam = rules.get_bam.output.bam,
        bed = config["HUMAN_TRANSCRIPT_BED_GZ"]
    output:
        bed1 = temp(OUTDIR + "/02_read_location/{sample}.bed"),
        bed2 = OUTDIR + "/02_read_location/{sample}.bed.gz",
        tsv = OUTDIR + "/02_read_location/{sample}_summary.tsv",
    log:
        OUTDIR + "/02_read_location/{sample}.log"
    threads:
        8
    shell:
        """
        stat_read_location.py -b {input.bam} -g {input.bed} -t {threads} -o {output.bed1} -s {output.tsv} &> {log}
        bgzip -c {output.bed1} > {output.bed2}
        """
