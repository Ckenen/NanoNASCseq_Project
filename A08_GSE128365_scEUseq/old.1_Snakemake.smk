#!/usr/bin/runsnakemake

# K562

SAMPLES = [
    "SRR11195416", # HeatShock-37C-NS-EU-45min-K562
    "SRR11195417", # HeatShock-37C-NS-EU-45min-K562
    "SRR11195418", # HeatShock-37C-Tot-EU-45min-K562
    "SRR11195419", # HeatShock-37C-Tot-EU-45min-K562
    "SRR11195420", # HeatShock-37C-Tot-EU-45min-K562
    "SRR11195421", # HeatShock-37C-Tot-EU-45min-K562
    "SRR11195422", # HeatShock-42C-Tot-DMSO-45min-K562
    "SRR11195423", # HeatShock-42C-Tot-DMSO-45min-K562
    "SRR11195424", # HeatShock-42C-Tot-DMSO-45min-K562
    "SRR11195425", # HeatShock-42C-Tot-DMSO-45min-K562
    "SRR11195426", # HeatShock-42C-NS-EU-45min-K562
    "SRR11195427", # HeatShock-42C-NS-EU-45min-K562
    "SRR11195428", # HeatShock-42C-NS-EU-45min-K562
    "SRR11195429", # HeatShock-42C-NS-EU-45min-K562
    "SRR11195430", # HeatShock-42C-Tot-EU-45min-K562
    "SRR11195431", # HeatShock-42C-Tot-EU-45min-K562
    "SRR11195432", # HeatShock-42C-Tot-EU-45min-K562
    "SRR11195433", # HeatShock-42C-Tot-EU-45min-K562
    "SRR11195434", # HeatShock-37C-Tot-DMSO-45min-K562
    "SRR11195435", # HeatShock-37C-Tot-DMSO-45min-K562
    "SRR11195436", # HeatShock-37C-Tot-DMSO-45min-K562
    "SRR11195437", # HeatShock-37C-Tot-DMSO-45min-K562
    "SRR11195438", # HeatShock-37C-NS-EU-45min-K562
    "SRR11195439", # HeatShock-37C-NS-EU-45min-K562
]
SAMPLES = SAMPLES[:5]


OUTDIR = "results"

rule all:
    input:
        expand(OUTDIR + "/sra/{sample}.sra", sample=SAMPLES),
        expand(OUTDIR + "/fastq/{sample}_1.fastq.gz", sample=SAMPLES),
        expand(OUTDIR + "/star/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/filtered/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/infer_experiment/{sample}.txt", sample=SAMPLES),

rule prefetch:
    output:
        sra = OUTDIR + "/sra/{sample}.sra"
    log:
        OUTDIR + "/sra/{sample}.log"
    conda:
        "sratools"
    shell:
        """
        prefetch --output-file {output.sra} {wildcards.sample} &> {log}
        """

rule sra2fq:
    input:
        sra = rules.prefetch.output.sra
    output:
        fq1 = temp(OUTDIR + "/fastq/{sample}_1.fastq"),
        fq2 = temp(OUTDIR + "/fastq/{sample}_2.fastq"),
        fq3 = OUTDIR + "/fastq/{sample}_1.fastq.gz",
        fq4 = OUTDIR + "/fastq/{sample}_2.fastq.gz",
    log:
        OUTDIR + "/fastq/{sample}.log"
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

rule star_mapping:
    input:
        fq = rules.sra2fq.output.fq4,
        idx = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.v39.STAR.index"
    output:
        out = directory(OUTDIR + "/star/{sample}")
    log:
        OUTDIR + "/star/{sample}.log"
    params:
        prefix = OUTDIR + "/star/{sample}/{sample}"
    threads:
        12
    conda:
        "star"
    shell:
        """(
        mkdir -p {output.out}
        STAR --runThreadN {threads} \
            --outFileNamePrefix {params.prefix}. \
            --genomeDir {input.idx} \
            --genomeLoad LoadAndKeep \
            --readFilesCommand zcat \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fq} ) &> {log}
        """

rule filter_bam:
    input:
        bamdir = rules.star_mapping.output.out
    output:
        bam = OUTDIR + "/filtered/{sample}.bam"
    log:
        OUTDIR + "/filtered/{sample}.log"
    shell:
        """(
        samtools view -@ {threads} -q 30 -d "NH:1" -F 2308 \
            --expr 'rname =~ "^chr([0-9]+|[XY])$"' -o {output.bam} \
            {input.bamdir}/{wildcards.sample}.Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule infer_experiment:
    input:
        bam = rules.filter_bam.output.bam,
        bed = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed"
    output:
        txt = OUTDIR + "/infer_experiment/{sample}.txt"
    log:
        OUTDIR + "/infer_experiment/{sample}.log"
    shell:
        """
        infer_experiment.py -s 2000000 -i {input.bam} -r {input.bed} > {output.txt} 2> {log}
        """