#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/rmdup"
OUTDIR = "results/snps"

rule all:
    input:
        # expand(OUTDIR + "/vcfs/{sample}.vcf.gz", sample=SAMPLES[:1]),
        # expand(OUTDIR + "/haplotag/{sample}.bam", sample=SAMPLES),
        # expand(OUTDIR + "/counts/{sample}.tsv", sample=SAMPLES),
        # OUTDIR + "/merge_strategy/all_samples_merged.bam",
        # OUTDIR + "/merge_strategy/all_samples_merged.vcf.gz",

rule call_het_snps:
    input:
        bam = INDIR + "/{sample}.human.bam"
    output:
        vcf = temp(OUTDIR + "/vcfs/{sample}.vcf"),
        vcf_gz = OUTDIR + "/vcfs/{sample}.vcf.gz"
    log:
        OUTDIR + "/vcfs/{sample}.log"
    threads:
        THREADS
    shell:
        """(
        ./scripts/call_het_snps_from_rnaseq.py \
            {input.bam} {wildcards.sample} {threads} {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz} ) &> {log}
        """

rule haplotag:
    input:
        vcf = OUTDIR + "/vcfs/20221128_K562_Actd_0h_rep1.vcf.gz",
        bam = INDIR + "/{sample}.human.bam",
        fsa = get_fasta("human")
    output:
        bam = OUTDIR + "/haplotag/{sample}.bam"
    conda:
        "whatshap"
    log:
        OUTDIR + "/haplotag/{sample}.log"
    shell:
        """
        whatshap haplotag --ignore-read-groups \
            -r {input.fsa} \
            -o {output.bam} \
            {input.vcf} {input.bam} &> {log}
        samtools index {output.bam}
        """

rule stat_hp_reads:
    input:
        vcf = OUTDIR + "/vcfs/20221128_K562_Actd_0h_rep1.vcf.gz",
        bam = OUTDIR + "/haplotag/{sample}.bam"
    output:
        tsv = OUTDIR + "/counts/{sample}.tsv"
    log:
        OUTDIR + "/counts/{sample}.log"
    shell:
        """
        ./scripts/stat_hp_reads.py {input.vcf} {input.bam} {output.tsv} &> {log}
        """

rule merge_bam:
    input:
        bams = expand(INDIR + "/{sample}.human.bam", sample=SAMPLES)
    output:
        bam = OUTDIR + "/merge_strategy/all_samples_merged.bam"
    log:
        OUTDIR + "/merge_strategy/all_samples_merged.log"
    threads:
        4
    shell:
        """(
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule call_het_snps_from_merged_bam:
    input:
        bam = OUTDIR + "/merge_strategy/all_samples_merged.bam"
    output:
        vcf = temp(OUTDIR + "/merge_strategy/all_samples_merged.vcf"),
        vcf_gz = OUTDIR + "/merge_strategy/all_samples_merged.vcf.gz"
    log:
        OUTDIR + "/merge_strategy/all_samples_merged.log"
    threads:
        THREADS
    shell:
        """(
        ./scripts/call_het_snps_from_rnaseq.py \
            {input.bam} AllSample {threads} {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz} ) &> {log}
        """