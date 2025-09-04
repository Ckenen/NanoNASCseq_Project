#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
BOWTIE2_INDEX = config["RIBO_BOWTIE2_INDEX"]
FQDIR = "data/datasets"
OUTDIR = "results/1_prepare"

rule all:
    input:
        # expand(OUTDIR + "/1_cutadapt/{sample}_R1.fastq.gz", sample=SAMPLES),
        expand(OUTDIR + "/2_bowtie2/{sample}.bam", sample=SAMPLES)

rule cutadapt:
    input:
        fq1 = FQDIR + "/{sample}_R1.fastq.gz",
        fq2 = FQDIR + "/{sample}_R2.fastq.gz"
    output:
        fq1 = OUTDIR + "/1_cutadapt/{sample}_R1.fastq.gz",
        fq2 = OUTDIR + "/1_cutadapt/{sample}_R2.fastq.gz"
    conda:
        "cutadapt"
    log:
        OUTDIR + "/1_cutadapt/{sample}.log"
    threads:
        12
    shell:
        """
        cutadapt --max-n 2 -m 20 -q 30 -j {threads} \
            -a AGATCGGAAGAGCACACGTC -a GATCGGAAGAGCACACGTCT \
            -A AGATCGGAAGAGCGTCGTGT -A GATCGGAAGAGCGTCGTGTA \
            -o {output.fq1} -p {output.fq2} \
            {input.fq1} {input.fq2} &> {log}
        """

rule bowtie2:
    input:
        fq1 = rules.cutadapt.output.fq1,
        fq2 = rules.cutadapt.output.fq2,
        idx = BOWTIE2_INDEX,
    output:
        bam = OUTDIR + "/2_bowtie2/{sample}.bam",
        fq1 = OUTDIR + "/2_bowtie2/{sample}.unmapped.fastq.1.gz",
        fq2 = OUTDIR + "/2_bowtie2/{sample}.unmapped.fastq.2.gz"
    log:
        OUTDIR + "/2_bowtie2/{sample}.log"
    conda:
        "bowtie2"
    params:
        prefix = OUTDIR + "/2_bowtie2/{sample}"
    threads:
        12
    shell:
        """(
        bowtie2 -p {threads} --local --no-unal --un-conc-gz {params.prefix}.unmapped.fastq.gz \
            -x {input.idx}/ref -1 {input.fq1} -2 {input.fq2} \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {params.prefix}_TMP - > {output.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """