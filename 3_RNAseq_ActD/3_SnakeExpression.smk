#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
SPECIES = config["SPECIES"]
BAMDIR = "results/2_mapping/4_rmdup"
OUTDIR = "results/3_expression"

rule all:
    input:
        expand(OUTDIR + "/1_fpkm/{sample}.{species}.tsv", sample=SAMPLES, species=SPECIES),
        expand(OUTDIR + "/2_intron_counts/{sample}.tsv", sample=SAMPLES),

rule calculate_fpkm:
    input:
        bam = BAMDIR + "/{sample}.{species}.bam",
        bed = lambda wildcards: config["HUMAN_TRANSCRIPT_BED"] if wildcards.species == "human" else config["FLY_TRANSCRIPT_BED"],
        tsv = lambda wildcards: config["HUMAN_ANNOTATION_CSV"] if wildcards.species == "human" else config["FLY_ANNOTATION_CSV"],
    output:
        tsv = OUTDIR + "/1_fpkm/{sample}.{species}.tsv"
    log:
        OUTDIR + "/1_fpkm/{sample}.{species}.log"
    threads:
        8
    shell:
        """
        nasctools CalculateFPKM -t {threads} -s R -l PE -a {input.tsv} {input.bam} {input.bed} {output.tsv} &> {log}
        """

rule quantify_introns:
    input:
        bam = BAMDIR + "/{sample}.human.bam",
        fa = config["HUMAN_FASTA"],
        bed = config["HUMAN_TRANSCRIPT_BED"]
    output:
        tsv = OUTDIR + "/2_intron_counts/{sample}.tsv"
    log:
        OUTDIR + "/2_intron_counts/{sample}.log"
    threads:
        8
    shell:
        """
        ./scripts/quantify_introns.py -t {threads} -d R -f {input.fa} -b {input.bed} {input.bam} {output.tsv} &> {log}
        """