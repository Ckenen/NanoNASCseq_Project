#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
OUTDIR = "results/01_mapping"

rule all:
    input:
        expand(OUTDIR + "/01_minimap2_mapped/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/02_bam_filtered/{sample}.bam", sample=SAMPLES),

rule minimap2:
    input:
        fq = "data/{sample}.fastq.gz",
        mmi = config["MINIMAP2_GENOME"],
        bed = config["TRANSCRIPT_BED"]
    output:
        bam = OUTDIR + "/01_minimap2_mapped/{sample}.bam"
    log:
        OUTDIR + "/01_minimap2_mapped/{sample}.log"
    params:
        rg = '@RG\\tID:{sample}\\tLB:{sample}\\tSM:{sample}'
    threads:
        24
    shell:
        """(
        minimap2 -ax splice -u f -Y --MD -R '{params.rg}' -t {threads} \
            --junc-bed {input.bed} {input.mmi} {input.fq} \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {output.bam} -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = OUTDIR + "/02_bam_filtered/{sample}.bam"
    log:
        OUTDIR + "/02_bam_filtered/{sample}.log"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            -q 30 -m 200 -F 2308 -o {output.bam} {input.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule sqanti3:
    input:
        gtf1 = rules.collapse.output.gtf,
        gtf2 = lambda wildcards: get_annotation_gtf(wildcards.cell),
        fa = lambda wildcards: get_genome_fasta(wildcards.cell)
    output:
        out = directory(OUTDIR + "/03_sqanti3/{sample}")
    log:
        OUTDIR + "/03_sqanti3/{sample}.log"
    conda:
        "SQANTI3.env"
    threads:
        4
    shell:
        """
        ./scripts/assembly/run_sqanti3_clean.sh {input} {threads} {output} &> {log}
        """