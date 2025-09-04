#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
SAMPLES = config["SAMPLES"]
SCRIPT_ROOT = config["SCRIPT_ROOT"]
INDIR = "results/02_dropseq/12_DetectedErrors"
OUTDIR = "results/03_mismatch"

rule all:
    input:
        expand(OUTDIR + "/01_BamToSam/{sample}.sam", sample=SAMPLES),
        expand(OUTDIR + "/02_AnnoRead/{sample}.txt", sample=SAMPLES),
        expand(OUTDIR + "/03_SamToTsv/{sample}.tsv", sample=SAMPLES),
        expand(OUTDIR + "/03_SamToTsv/{sample}.tsv_q27_gene_anno_stat.txt", sample=SAMPLES),
        expand(OUTDIR + "/04_TC/{sample}.tsv", sample=SAMPLES),
        expand(OUTDIR + "/04_TC/{sample}.tsv_q27.tsv", sample=SAMPLES),

rule BamToSam:
    input:
        bam = INDIR + "/{sample}.bam"
    output:
        sam = OUTDIR + "/01_BamToSam/{sample}.sam"
    threads:
        8
    shell:
        """
        samtools view -@ {threads} -h {input.bam} > {output.sam}
        """

rule GetAnnoRead:
    input:
        sam = rules.BamToSam.output.sam
    output:
        txt = OUTDIR + "/02_AnnoRead/{sample}.txt"
    shell:
        """
        grep "GE:Z" {input.sam} | awk '{{print $1}}' > {output.txt}
        """

rule SamToTsv:
    input:
        sam = rules.BamToSam.output.sam,
        fa = lambda wildcards: get_genome_fasta(wildcards.sample)
    output:
        tsv = OUTDIR + "/03_SamToTsv/{sample}.tsv"
    log:
        OUTDIR + "/03_SamToTsv/{sample}.log"
    shell:
        """
        ( sam2tsv --reference {input.fa} {input.sam} | awk '{{if ($9 ~/M|=|X/ ) print $0}}' > {output.tsv} ) &> {log}
        """

rule extrac_conversion_frequency_gene_annotate:
    input:
        tsv1 = rules.GetAnnoRead.output.txt,
        tsv2 = rules.SamToTsv.output.tsv
    output:
        tsv = OUTDIR + "/03_SamToTsv/{sample}.tsv_q27_gene_anno_stat.txt"
    shell:
        """
        perl {SCRIPT_ROOT}/extrac_conversion_frequency_gene_annotate.pl \
            -read {input.tsv1} \
            -tsv {input.tsv2} \
            -qual 27
        """

rule GetTC:
    input:
        tsv = rules.SamToTsv.output.tsv
    output:
        tsv = OUTDIR + "/04_TC/{sample}.tsv"
    shell:
        """
        awk '{{if ($2 ==0 && $5=="C" && $8=="T") print $0; else if ( $2==16 && $5=="G" && $8=="A" ) print $0}}' {input.tsv} > {output.tsv}
        """

rule FilterTC:
    input:
        tsv = rules.GetTC.output.tsv
    output:
        tsv = OUTDIR + "/04_TC/{sample}.tsv_q27.tsv"
    shell:
        """
        perl {SCRIPT_ROOT}/extrac_refT_readC.pl -tsv {input.tsv} -qual 27
        """

