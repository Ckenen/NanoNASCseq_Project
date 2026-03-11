#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
SPECIES = ["human", "fly"]
GROUPS = ["filtered", "rmdup"]
FQDIR = "data/datasets"
OUTDIR = "results/1_rnaseq"

rule all:
    input:
        OUTDIR + "/0_index/ribosomal_bowtie2.idx",
        OUTDIR + "/0_index/genome_star.idx",
        # expand(OUTDIR + "/1_cutadapt/{sample}_R1.fastq.gz", sample=SAMPLES),
        expand(OUTDIR + "/2_bowtie2/{sample}.unmapped.fastq.1.gz", sample=SAMPLES),
        expand(OUTDIR + "/3_star/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand(OUTDIR + "/4_bams/{sample}.{species}.{group}.bam", sample=SAMPLES, species=SPECIES, group=GROUPS),
        expand(OUTDIR + "/4_bams/{sample}.{species}.{group}.flagstat", sample=SAMPLES, species=SPECIES, group=GROUPS),
        expand(OUTDIR + "/4_bams/{sample}.{species}.{group}.fpkm.tsv", sample=SAMPLES, species=SPECIES, group=GROUPS),
        expand(OUTDIR + "/4_bams/{sample}.{species}.{group}.read_location.tsv", sample=SAMPLES, species=SPECIES, group=GROUPS),
        expand(OUTDIR + "/4_bams/{sample}.{species}.{group}.intron_count.tsv", sample=SAMPLES, species=SPECIES, group=GROUPS),

# Building index

rule bowtie2_build:
    input:
        fa = config["FASTA_RIBO"]
    output:
        idx = directory(OUTDIR + "/0_index/ribosomal_bowtie2.idx")
    log:
        OUTDIR + "/0_index/ribosomal_bowtie2.idx.log"
    conda:
        "bowtie2"
    shell:
        """(
        mkdir {output.idx}
        bowtie2-build {input.fa} {output.idx}/ref ) &> {log}
        """

rule star_build:
    input:
        fa1 = config["HUMAN_FASTA"],
        fa2 = config["FLY_FASTA"],
        gtf1 = config["HUMAN_GTF"],
        gtf2 = config["FLY_GTF"]
    output:
        fa = OUTDIR + "/0_index/genome_star.fa",
        fai = OUTDIR + "/0_index/genome_star.fa.fai",
        gtf = OUTDIR + "/0_index/genome_star.gtf",
        gtf2 = OUTDIR + "/0_index/genome_star.gtf.gz",
        tbi2 = OUTDIR + "/0_index/genome_star.gtf.gz.tbi",
        idx = directory(OUTDIR + "/0_index/genome_star.idx")
    log:
        OUTDIR + "/0_index/genome_star.log"
    conda:
        "star"
    threads:
        24
    shell:
        """
        cat {input.fa1} {input.fa2} > {output.fa}
        samtools faidx {output.fa}
        cat {input.gtf1} {input.gtf2} | grep -v '#' | sort -k1,1 -k4,4n -k5,5n > {output.gtf} 
        bgzip -c {output.gtf} > {output.gtf2}
        tabix -p gff {output.gtf2}
        mkdir -p {output.idx}
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.idx} \
            --genomeFastaFiles {output.fa} --sjdbGTFfile {output.gtf} &> {log}
        """

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
        idx = rules.bowtie2_build.output.idx,
    output:
        bam = temp(OUTDIR + "/2_bowtie2/{sample}.bam"),
        bai = temp(OUTDIR + "/2_bowtie2/{sample}.bam.bai"),
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

rule star:
    input:
        fq1 = rules.bowtie2.output.fq1,
        fq2 = rules.bowtie2.output.fq2,
        idx = rules.star_build.output.idx
    output:
        bam = OUTDIR + "/3_star/{sample}.Aligned.sortedByCoord.out.bam"
    log:
        OUTDIR + "/3_star/{sample}.log"
    params:
        prefix = OUTDIR + "/3_star/{sample}."
    conda:
        "star"
    threads:
        12
    shell:
        """(
        mkdir -p {output}
        STAR --runThreadN {threads} \
            --outFileNamePrefix {params.prefix} \
            --genomeDir {input.idx} \
            --readFilesCommand zcat \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fq1} {input.fq2} ) &> {log}
        """

rule filter_and_split: # filter and split
    input:
        bam = rules.star.output.bam
    output:
        bam = OUTDIR + "/4_bams/{sample}.{species}.filtered.bam",
        bai = OUTDIR + "/4_bams/{sample}.{species}.filtered.bam.bai"
    log:
        OUTDIR + "/4_bams/{sample}.{species}.filtered.log"
    params:
        pattern = lambda wildcards: config["%s_SEQNAME_PATTERN" % wildcards.species.upper()]
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} -q 30 -d "NH:1" -f 2 -F 2308 \
            --expr 'rname =~ "{params.pattern}"' -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule rmdup:
    input:
        bam = rules.filter_and_split.output.bam
    output:
        bam = OUTDIR + "/4_bams/{sample}.{species}.rmdup.bam",
        bai = OUTDIR + "/4_bams/{sample}.{species}.rmdup.bam.bai"
    log:
        OUTDIR + "/4_bams/{sample}.{species}.rmdup.log"
    threads:
        4
    shell:
        """(
        sambamba markdup -r -t {threads} {input.bam} {output.bam} 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule calculate_fpkm:
    input:
        bam = OUTDIR + "/4_bams/{sample}.{species}.{group}.bam",
        bed = lambda wildcards: config["%s_TRANSCRIPT_BED" % wildcards.species.upper()],
        tsv = lambda wildcards: config["%s_ANNOTATION_CSV" % wildcards.species.upper()]
    output:
        tsv = OUTDIR + "/4_bams/{sample}.{species}.{group}.fpkm.tsv"
    log:
        OUTDIR + "/4_bams/{sample}.{species}.{group}.fpkm.log"
    threads:
        8
    shell:
        """
        nasctools CalculateFPKM -t {threads} -s R -a {input.tsv} {input.bam} {input.bed} {output.tsv} &> {log}
        """

rule quantify_introns:
    input:
        bam = OUTDIR + "/4_bams/{sample}.{species}.{group}.bam",
        fa = lambda wildcards: config["%s_FASTA" % wildcards.species.upper()],
        bed = lambda wildcards: config["%s_TRANSCRIPT_BED" % wildcards.species.upper()]
    output:
        tsv = OUTDIR + "/4_bams/{sample}.{species}.{group}.intron_count.tsv"
    log:
        OUTDIR + "/4_bams/{sample}.{species}.{group}.intron_count.log"
    threads:
        8
    shell:
        """
        ./scripts/quantify_introns.py -t {threads} -d R -f {input.fa} -b {input.bed} {input.bam} {output.tsv} &> {log}
        """

rule stat_read_location:
    input:
        bam = OUTDIR + "/4_bams/{sample}.{species}.{group}.bam",
        bed = lambda wildcards: config["%s_TRANSCRIPT_BED" % wildcards.species.upper()]
    output:
        bed1 = temp(OUTDIR + "/4_bams/{sample}.{species}.{group}.read_location.bed"),
        bed2 = OUTDIR + "/4_bams/{sample}.{species}.{group}.read_location.bed.gz",
        tsv = OUTDIR + "/4_bams/{sample}.{species}.{group}.read_location.tsv"
    log:
        OUTDIR + "/4_bams/{sample}.{species}.{group}.read_location.log"
    threads:
        8
    shell:
        """
        stat_read_location.py -b {input.bam} -g {input.bed} -t {threads} -o {output.bed1} -s {output.tsv} &> {log}
        bgzip -c {output.bed1} > {output.bed2}
        """

# Common rules

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """