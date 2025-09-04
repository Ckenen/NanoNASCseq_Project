#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
FQDIR = "results/01_prepare/02_bowtie2"
OUTDIR = "results/02_mapping"

rule all:
    input:
        # expand(OUTDIR + "/01_star/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/02_filtered/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/02_filtered/{sample}.flagstat", sample=SAMPLES),
        expand(OUTDIR + "/03_rmdup/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/03_rmdup/{sample}.flagstat", sample=SAMPLES),

rule star_mapping:
    input:
        fq = FQDIR + "/{sample}.unmapped.fastq.gz",
        idx = config["GENOME_STAR_INDEX"]
    output:
        out = directory(OUTDIR + "/01_star/{sample}.out"),
        bam = temp(OUTDIR + "/01_star/{sample}.bam"),
        bai = temp(OUTDIR + "/01_star/{sample}.bam.bai")
    conda:
        "star"
    log:
        OUTDIR + "/01_star/{sample}.log"
    threads:
        12
    shell:
        """(
        mkdir -p {output.out}
        STAR --runThreadN {threads} \
            --outFileNamePrefix {output.out}/{wildcards.sample}. \
            --genomeDir {input.idx} \
            --readFilesCommand zcat \
            --alignEndsType EndToEnd \
            --outFilterMismatchNoverLmax 0.05 \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fq}
        mv {output.out}/{wildcards.sample}.Aligned.sortedByCoord.out.bam {output.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.star_mapping.output.bam
    output:
        bam = OUTDIR + "/02_filtered/{sample}.bam",
        bai = OUTDIR + "/02_filtered/{sample}.bam.bai"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} -q 30 -d "NH:1" --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            -F 2308 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule rmdup:
    input:
        bam = rules.filter_bam.output.bam
    output:
        tmpdir = temp(directory(OUTDIR + "/03_rmdup/{sample}.TMP")),
        bam = OUTDIR + "/03_rmdup/{sample}.bam",
        bai = OUTDIR + "/03_rmdup/{sample}.bam.bai"
    log:
        OUTDIR + "/03_rmdup/{sample}.log"
    threads:
        4
    shell:
        """
        mkdir -p {output.tmpdir}
        sambamba markdup -r -t {threads} --tmpdir={output.tmpdir} {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

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
        