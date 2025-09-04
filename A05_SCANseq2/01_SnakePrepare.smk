#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
SRRS = [sample.split("_")[0] for sample in SAMPLES]
OUTDIR = "results/01_prepare"

rule all:
    input:
        # expand(OUTDIR + "/01_sra/{srr}.sra", srr=SRRS),
        expand(OUTDIR + "/02_fastq/{sample}.fastq.gz", sample=SAMPLES),
        expand(OUTDIR + "/03_read_count/{sample}.txt", sample=SAMPLES),
        
rule prefetch:
    output:
        sra = OUTDIR + "/01_sra/{srr}.sra"
    log:
        OUTDIR + "/01_sra/{srr}.log"
    conda:
        "sratools"
    shell:
        """
        prefetch --max-size 200000000 -o {output} {wildcards.srr} &> {log}
        """

rule sra2fastq:
    input:
        sra = lambda wildcards: OUTDIR + "/01_sra/%s.sra" % wildcards.sample.split("_")[0]
    output:
        fq = OUTDIR + "/02_fastq/{sample}.fastq.gz"
    log:
        OUTDIR + "/02_fastq/{sample}.log"
    threads:
        12
    conda:
        "sratools"
    shell:
        """(
        fasterq-dump -e {threads} -O `dirname {output.fq}` {input.sra}
        pigz -p {threads} `dirname {output.fq}`/`basename {output.fq} .gz` ) &> {log}
        """

rule get_read_count:
    input:
        fq = rules.sra2fastq.output.fq
    output:
        txt = OUTDIR + "/03_read_count/{sample}.txt"
    shell:
        """
        zcat {input.fq} | wc -l | awk '{{print $1/4}}' > {output.txt}
        """