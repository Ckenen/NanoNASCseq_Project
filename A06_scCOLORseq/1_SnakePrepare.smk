#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
SRRS = [s.split("_")[0] for s in SAMPLES]
OUTDIR = "results/01_prepare"

rule all:
    input:
        expand(OUTDIR + "/01_sra/{srr}.sra", srr=SRRS),
        expand(OUTDIR + "/02_fastq/{sample}.fastq.gz", sample=SAMPLES),
        expand(OUTDIR + "/03_split/{sample}", sample=SAMPLES),

rule prefetch:
    output:
        sra = OUTDIR + "/01_sra/{srr}.sra"
    log:
        OUTDIR + "/01_sra/{srr}.log"
    conda:
        "sratools"
    shell:
        """
        prefetch --max-size 200000000 -o {output.sra} {wildcards.srr} &> {log}
        """
    
rule sra2fastq:
    input:
        sra = lambda wildcards: OUTDIR + "/01_sra/%s.sra" % wildcards.sample.split("_")[0]
    output:
        fq1 = temp(OUTDIR + "/02_fastq/{sample}.fastq"),
        fq2 = OUTDIR + "/02_fastq/{sample}.fastq.gz"
    log:
        OUTDIR + "/02_fastq/{sample}.log"
    threads:
        6
    conda:
        "sratools"
    shell:
        """(
        fasterq-dump -e {threads} -o {output.fq1} {input.sra}
        pigz -p {threads} -c {output.fq1} > {output.fq2} ) &> {log}
        """

rule split:
    input:
        fq = rules.sra2fastq.output.fq2
    output:
        directory(OUTDIR + "/03_split/{sample}")
    log:
        OUTDIR + "/03_split/{sample}.log"
    threads:
        8
    shell:
        """(
        mkdir -p {output}
        zcat {input.fq} | split --additional-suffix=.fastq -l 1000000 - {output}/
        pigz -p {threads} {output}/*.fastq ) &> {log}
        """