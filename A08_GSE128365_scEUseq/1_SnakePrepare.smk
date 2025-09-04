#!/usr/bin/env runsnakemake
import pandas as pd
DAT = pd.read_csv("data/SraRunTable.PRJNA527268.csv")
DAT = DAT[DAT["cell_type"] == "RPE1-FUCCI"]
RUNS = list(DAT["Run"])
OUTDIR = "results/01_prepare"

rule all:
    input:
        expand(OUTDIR + "/01_sra/{run}.sra", run=RUNS),
        expand(OUTDIR + "/02_fastq/{run}_1.fastq.gz", run=RUNS),

rule prefetch:
    output:
        sra = OUTDIR + "/01_sra/{run}.sra"
    log:
        OUTDIR + "/01_sra/{run}.log"
    conda:
        "sratools"
    shell:
        """
        prefetch --output-file {output.sra} {wildcards.run} &> {log}
        """

rule sra2fq:
    input:
        sra = rules.prefetch.output.sra
    output:
        fq1 = temp(OUTDIR + "/02_fastq/{run}_1.fastq"),
        fq2 = temp(OUTDIR + "/02_fastq/{run}_2.fastq"),
        fq3 = OUTDIR + "/02_fastq/{run}_1.fastq.gz",
        fq4 = OUTDIR + "/02_fastq/{run}_2.fastq.gz",
    log:
        OUTDIR + "/02_fastq/{run}.log"
    threads:
        6
    conda:
        "sratools"
    shell:
        """(
        fasterq-dump {input.sra} --outdir `dirname {output.fq1}` --threads {threads}
        pigz -p {threads} -c {output.fq1} > {output.fq3}
        pigz -p {threads} -c {output.fq2} > {output.fq4} ) &> {log}
        """