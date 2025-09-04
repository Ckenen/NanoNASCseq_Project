#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"][:5]
# SAMPLES = ["GSM3672480_SRR8733760_RPE1_Chase_plate01_unlabeled"]
FQDIR = "results/01_prepare/02_fastq"
OUTDIR = "results/02_mapping"
GENOME_FASTA = "../1_NanoNASCseq/data/GRCh38_FUCCI.fa"
ANNOTATION_GTF = "../1_NanoNASCseq/data/gencode.v39.annotation.fucci.gtf"

rule all:
    input:
        expand(OUTDIR + "/01_tag_bc_umi/{sample}.fastq.gz", sample=SAMPLES),
        expand(OUTDIR + "/02_cutadapt/{sample}.fastq.gz", sample=SAMPLES),
        OUTDIR + "/03_index/star_index",
        expand(OUTDIR + "/04_star/{sample}", sample=SAMPLES),

rule tag_bc_umi:
    input:
        fq1 = lambda wildcards: FQDIR + "/%s_1.fastq.gz" % wildcards.sample.split("_")[1],
        fq2 = lambda wildcards: FQDIR + "/%s_2.fastq.gz" % wildcards.sample.split("_")[1],
    output:
        fq = OUTDIR + "/01_tag_bc_umi/{sample}.fastq.gz"
    log:
        OUTDIR + "/01_tag_bc_umi/{sample}.log"
    shell:
        """
        ./scripts/tag_bc_umi.py {input.fq1} {input.fq2} 2> {log} | gzip -c > {output.fq}
        """

rule cutadapt:
    input:
        fq = rules.tag_bc_umi.output.fq
    output:
        fq = OUTDIR + "/02_cutadapt/{sample}.fastq.gz"
    log:
        OUTDIR + "/02_cutadapt/{sample}.log"
    threads:
        12
    conda:
        "cutadapt"
    shell:
        """
        cutadapt -j {threads} -m 20 -a AAAAAAAAAAAAAAAAAAAA -o {output.fq} {input.fq} &> {log}
        """

rule star_index:
    input:
        fa = GENOME_FASTA,
        gtf = ANNOTATION_GTF
    output:
        directory(OUTDIR + "/03_index/star_index")
    log:
        OUTDIR + "/03_index/star_index.log"
    threads:
        20
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} &> {log}
        """

rule star:
    input:
        fq = rules.cutadapt.output.fq,
        idx = rules.star_index.output
    output:
        directory(OUTDIR + "/04_star/{sample}")
    log:
        OUTDIR + "/04_star/{sample}.log"
    threads:
        20
    shell:
        """
        mkdir -p {output}
        STAR --runThreadN {threads} \
            --genomeDir {input.idx} \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesCommand zcat \
            --readFilesIn {input.fq} \
            --outFileNamePrefix {output}/{wildcards.sample}. &> {log}
        samtools index -@ 4 {output}/{wildcards.sample}.Aligned.sortedByCoord.out.bam
        """