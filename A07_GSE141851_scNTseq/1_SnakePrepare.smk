#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
SRAS = [s.split("_")[0] for s in SAMPLES]
OUTDIR = "results/01_prepare"

rule all:
    input:
        # expand(OUTDIR + "/01_sra/{sra}.sra", sra=SRAS),
        expand(OUTDIR + "/02_fastq/{sample}_1.fastq.gz", sample=SAMPLES),

rule prefetch:
    output:
        sra = temp(OUTDIR + "/01_sra/{sra}.sra")
    log:
        log = OUTDIR + "/01_sra/{sra}.log"
    conda:
        "sratools"
    shell:
        """
        prefetch -o {output.sra} {wildcards.sra} &> {log}
        """

rule sra2fq:
    input:
        sra = lambda wildcards: OUTDIR + "/01_sra/%s.sra" % wildcards.sample.split("_")[0]
    output:
        fq1 = OUTDIR + "/02_fastq/{sample}_1.fastq.gz",
        fq2 = OUTDIR + "/02_fastq/{sample}_2.fastq.gz"
    log:
        OUTDIR + "/02_fastq/{sample}.log"
    params:
        fq1 = lambda wildcards: OUTDIR + "/02_fastq/%s_1.fastq" % wildcards.sample.split("_")[0],
        fq2 = lambda wildcards: OUTDIR + "/02_fastq/%s_2.fastq" % wildcards.sample.split("_")[0]
    threads:
        6
    conda:
        "sratools"
    shell:
        """(
        fasterq-dump -3 -e {threads} -O `dirname {output.fq1}` {input.sra}
        pigz -p {threads} -c {params.fq1} > {output.fq1}
        pigz -p {threads} -c {params.fq2} > {output.fq2}
        rm {params.fq1} {params.fq2} ) &> {log}
        """