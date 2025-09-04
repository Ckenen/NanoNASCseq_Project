#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config['SAMPLES']
FQDIR = "data"
OUTDIR = "results/01_prepare"

rule all:
    input:
        expand(OUTDIR + "/01_cutadapt/{sample}_R1.fastq.gz", sample=SAMPLES),
        OUTDIR + "/02_bowtie2_index_ribo/ribo_bt2_index",
        expand(OUTDIR + "/03_bowtie2_mapping_ribo/{sample}.unmapped.fastq.1.gz", sample=SAMPLES),

rule cutadapt:
    input:
        fq1 = FQDIR + "/{sample}_R1.fastq.gz",
        fq2 = FQDIR + "/{sample}_R2.fastq.gz"
    output:
        fq1 = OUTDIR + "/01_cutadapt/{sample}_R1.fastq.gz",
        fq2 = OUTDIR + "/01_cutadapt/{sample}_R2.fastq.gz"
    log:
        OUTDIR + "/01_cutadapt/{sample}.log"
    threads:
        12
    conda:
        "cutadapt"
    shell:
        """
        cutadapt -m 20 -q 30 -j {threads} \
            -a AGTCTTGGAAGCTTG -A CTGGCTGCGACATCT \
            -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} &> {log}
        """

rule bowtie2_build_ribo:
    input:
        fa = config["RIBO_FASTA"]
    output:
        idx = directory(OUTDIR + "/02_bowtie2_index_ribo/ribo_bt2_index")
    log:
        OUTDIR + "/02_bowtie2_index_ribo/ribo_bt2_index.log"
    conda:
        "bowtie2"
    shell:
        """
        mkdir -p {output.idx}
        bowtie2-build --threads {threads} {input.fa} {output.idx}/ref &> {log}
        """

rule bowtie2_mapping_ribo:
    input:
        fq1 = rules.cutadapt.output.fq1,
        fq2 = rules.cutadapt.output.fq2,
        idx = rules.bowtie2_build_ribo.output.idx
    output:
        bam = OUTDIR + "/03_bowtie2_mapping_ribo/{sample}.bam",
        fq1 = OUTDIR + "/03_bowtie2_mapping_ribo/{sample}.unmapped.fastq.1.gz",
        fq2 = OUTDIR + "/03_bowtie2_mapping_ribo/{sample}.unmapped.fastq.2.gz"
    log:
        OUTDIR + "/03_bowtie2_mapping_ribo/{sample}.log"
    conda:
        "bowtie2"
    params:
        prefix = OUTDIR + "/03_bowtie2_mapping_ribo/{sample}"
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
