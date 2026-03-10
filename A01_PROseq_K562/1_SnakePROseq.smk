#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
GROUPS = ["allread", "rmdup"]
OUTDIR = "results/1_proseq"

rule all:
    input:
        # OUTDIR + "/0_index/ribosomal_bowtie2.idx",
        # OUTDIR + "/0_index/genome_star.idx",
        # expand(OUTDIR + "/1_cutadapt/{sample}.fastq.gz", sample=SAMPLES),
        # expand(OUTDIR + "/2_bowtie2/{sample}.unmapped.fastq.gz", sample=SAMPLES),
        # expand(OUTDIR + "/3_star/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand(OUTDIR + "/4_bams/{sample}.{group}.bam", sample=SAMPLES, group=GROUPS),
        expand(OUTDIR + "/4_bams/{sample}.{group}.flagstat", sample=SAMPLES, group=GROUPS),
        expand(OUTDIR + "/5_bigwigs/{sample}.{group}_raw_both.bw", sample=SAMPLES, group=GROUPS),
        expand(OUTDIR + "/6_bigwigs_single_base/{sample}.{group}_raw_both.bw", sample=SAMPLES, group=GROUPS),
        expand(OUTDIR + "/7_stats/{sample}.{group}.read_location.tsv", sample=SAMPLES, group=GROUPS),

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
        fa = config["FASTA"],
        gtf = config["GTF"]
    output:
        idx = directory(OUTDIR + "/0_index/genome_star.idx")
    log:
        OUTDIR + "/0_index/genome_star.log"
    conda:
        "star"
    threads:
        24
    shell:
        """
        mkdir -p {output.idx}
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.idx} \
            --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} &> {log}
        """

# Pipeline

rule cutadapt:
    input:
        fq = config["FQDIR"] + "/{sample}.fastq.gz"
    output:
        fq = OUTDIR + "/1_cutadapt/{sample}.fastq.gz"
    log:
        OUTDIR + "/1_cutadapt/{sample}.log"
    threads:
        12
    conda:
        "cutadapt"
    shell:
        """
        cutadapt -m 20 -q 30 -j {threads} -a TGGAATTCTCGGGTGCCAAGG -o {output.fq} {input.fq} &> {log}
        """

rule bowtie2:
    input:
        fq = rules.cutadapt.output.fq,
        idx = rules.bowtie2_build.output.idx
    output:
        fq = OUTDIR + "/2_bowtie2/{sample}.unmapped.fastq.gz"
    log:
        OUTDIR + "/2_bowtie2/{sample}.log"
    conda:
        "bowtie2"
    params:
        prefix = OUTDIR + "/2_bowtie2/{sample}"
    threads:
        12
    shell:
        """
        bowtie2 -p {threads} --local --un-gz {output.fq} -x {input.idx}/ref -U {input.fq} > /dev/null 2> {log}
        """

rule star:
    input:
        fq = rules.bowtie2.output.fq,
        idx = rules.star_build.output.idx
    output:
        bam = OUTDIR + "/3_star/{sample}.Aligned.sortedByCoord.out.bam",
        bai = OUTDIR + "/3_star/{sample}.Aligned.sortedByCoord.out.bam.bai"
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
        STAR --runThreadN {threads} \
            --outFileNamePrefix {params.prefix} \
            --genomeDir {input.idx} \
            --readFilesCommand zcat \
            --alignEndsType EndToEnd \
            --outFilterMismatchNoverLmax 0.05 \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fq}
        samtools index -@ 4 {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.star.output.bam
    output:
        bam = OUTDIR + "/4_bams/{sample}.allread.bam",
        bai = OUTDIR + "/4_bams/{sample}.allread.bam.bai"
    log:
        OUTDIR + "/4_bams/{sample}.allread.log"
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} -q 30 -d "NH:1" --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            -F 2308 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule rmdup:
    input:
        bam = rules.filter_bam.output.bam
    output:
        tmpdir = temp(directory(OUTDIR + "/4_bams/{sample}.rmdup.TMP")),
        bam = OUTDIR + "/4_bams/{sample}.rmdup.bam",
        bai = OUTDIR + "/4_bams/{sample}.rmdup.bam.bai"
    log:
        OUTDIR + "/4_bams/{sample}.rmdup.log"
    threads:
        4
    shell:
        """(
        mkdir -p {output.tmpdir}
        sambamba markdup -r -t {threads} --tmpdir={output.tmpdir} {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule bam2bw:
    input:
        bam = OUTDIR + "/4_bams/{sample}.{group}.bam"
    output:
        bw = OUTDIR + "/5_bigwigs/{sample}.{group}_raw_both.bw"
    log:
        OUTDIR + "/5_bigwigs/{sample}.{group}_raw_both.log"
    params:
        prefix = OUTDIR + "/5_bigwigs/{sample}.{group}"
    shell:
        """
        ../share/scripts/bam2bw.sh {input.bam} {params.prefix}
        """

rule bam2bw_single_base:
    input:
        bam = OUTDIR + "/4_bams/{sample}.{group}.bam"
    output:
        bw = OUTDIR + "/6_bigwigs_single_base/{sample}.{group}_raw_both.bw"
    log:
        OUTDIR + "/6_bigwigs_single_base/{sample}.{group}_raw_both.log"
    params:
        prefix = OUTDIR + "/6_bigwigs_single_base/{sample}.{group}"
    shell:
        """
        ../share/scripts/bam2bw_proseq.sh {input.bam} {params.prefix}
        """

rule stat_read_location:
    input:
        bam = OUTDIR + "/4_bams/{sample}.{group}.bam",
        bed = config["TRANSCRIPT_BED_GZ"]
    output:
        bed1 = temp(OUTDIR + "/7_stats/{sample}.{group}.read_location.bed"),
        bed2 = OUTDIR + "/7_stats/{sample}.{group}.read_location.bed.gz",
        tsv = OUTDIR + "/7_stats/{sample}.{group}.read_location.tsv"
    log:
        OUTDIR + "/7_stats/{sample}.{group}.read_location.log"
    threads:
        8
    shell:
        """(
        stat_read_location.py -b {input.bam} -g {input.bed} -t {threads} -o {output.bed1} -s {output.tsv}
        bgzip -c {output.bed1} > {output.bed2} ) &> {log}
        """

# common rules

rule flagstat:
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


