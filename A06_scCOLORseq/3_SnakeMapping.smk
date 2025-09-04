#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
FQDIR = "results/02_tallynn/13_remove_polya"
OUTDIR = "results/03_mapping"

rule all:
    input:
        OUTDIR + "/01_minimap2_index/genome.splice.mmi",
        expand(OUTDIR + "/02_minimap2_mapped/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/02_minimap2_mapped/{sample}.flagstat", sample=SAMPLES),
        expand(OUTDIR + "/03_bam_filtered/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/03_bam_filtered/{sample}.flagstat", sample=SAMPLES),
        expand(OUTDIR + "/04_stat_clip/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/04_stat_clip/{sample}.flagstat", sample=SAMPLES),

rule minimap2_build_index:
    input:
        fa = config["GENOME_FASTA"]
    output:
        mmi = OUTDIR + "/01_minimap2_index/genome.splice.mmi"
    log:
        OUTDIR + "/01_minimap2_index/genome.splice.log"
    threads:
        12
    shell:
        """
        minimap2 -t {threads} -x splice -d {output.mmi} {input.fa} &> {log}
        """

rule minimap2:
    input:
        fq = FQDIR + "/{sample}.fastq.gz",
        mmi = rules.minimap2_build_index.output.mmi,
        bed = config["TRANSCRIPT_BED"]
    output:
        bam = OUTDIR + "/02_minimap2_mapped/{sample}.bam"
    log:
        OUTDIR + "/02_minimap2_mapped/{sample}.log"
    threads:
        24
    shell: 
        """(
        minimap2 -ax splice -u f -Y --MD --junc-bed {input.bed} -t {threads} {input.mmi} {input.fq} \
            | samtools view -q 60 -@ {threads} -u - \
            | samtools sort -@ {threads} -T {output.bam} -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = OUTDIR + "/03_bam_filtered/{sample}.bam"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} --expr 'rname =~ "^(hg|mm)_chr([0-9]+|[XY])$"' -q 30 -m 200 -F 2308 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule stat_clip:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/04_stat_clip/{sample}.bam",
        txt = OUTDIR + "/04_stat_clip/{sample}.tsv"
    log:
        OUTDIR + "/04_stat_clip/{sample}.log"
    threads:
        4
    shell:
        """
        stat_clip.py -c 5 -s {output.txt} -o {output.bam} {input.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

# Common rules

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    shell:
        """
        samtools flagstat {input} > {output}
        """