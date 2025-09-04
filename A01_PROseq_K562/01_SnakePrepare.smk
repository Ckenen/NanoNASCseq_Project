#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
FQDIR = "data"
OUTDIR = "results/01_prepare"

rule all:
    input:
        # expand(OUTDIR + "/01_cutadapt/{sample}.fastq.gz", sample=SAMPLES),
        expand(OUTDIR + "/02_bowtie2/{sample}.unmapped.fastq.gz", sample=SAMPLES),

rule cutadapt:
    input:
        fq = FQDIR + "/{sample}.fastq.gz"
    output:
        fq = temp(OUTDIR + "/01_cutadapt/{sample}.fastq.gz")
    log:
        OUTDIR + "/01_cutadapt/{sample}.log"
    threads:
        12
    conda:
        "cutadapt"
    shell:
        """
        cutadapt -m 20 -q 30 -j {threads} -a TGGAATTCTCGGGTGCCAAGG \
            -o {output.fq} {input.fq} &> {log}
        """

rule bowtie2:
    input:
        fq = rules.cutadapt.output.fq,
        idx = config["RIBO_BOWTIE2_INDEX"]
    output:
        bam = OUTDIR + "/02_bowtie2/{sample}.bam",
        fq = OUTDIR + "/02_bowtie2/{sample}.unmapped.fastq.gz"
    log:
        OUTDIR + "/02_bowtie2/{sample}.log"
    conda:
        "bowtie2"
    params:
        prefix = OUTDIR + "/02_bowtie2/{sample}"
    threads:
        12
    shell:
        """(
        bowtie2 -p {threads} --local --no-unal --un-gz {output.fq} \
            -x {input.idx}/ref -U {input.fq} \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {params.prefix}_TMP - > {output.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """
