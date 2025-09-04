#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
BAMDIR = "results/01_mapping/02_bam_filtered"
OUTDIR = "results/03_expression"

rule all:
    input:
        expand(OUTDIR + "/02_collapsed/{sample}.gtf", sample=SAMPLES),
        expand(OUTDIR + "/03_sqanti3/{sample}", sample=SAMPLES),

rule collapse:
    input:
        bam = BAMDIR + "/{sample}.bam"
    output:
        bed = OUTDIR + "/02_collapsed/{sample}.bed",
        bed_gz = OUTDIR + "/02_collapsed/{sample}.bed.gz",
        gtf = OUTDIR + "/02_collapsed/{sample}.gtf",
        gtf_gz = OUTDIR + "/02_collapsed/{sample}.gtf.gz"
    log:
        OUTDIR + "/02_collapsed/{sample}.log"
    shell:
        """(
        ./scripts/expression/collapse_umi.py {input.bam} | sort -k1,1 -k2,2n -k3,3n > {output.bed}
        bgzip -c {output.bed} > {output.bed_gz}
        tabix -p bed -f {output.bed_gz}
        ./scripts/expression/bed2gtf.py {output.bed} | sort -k1,1 -k4,4n > {output.gtf}
        bgzip -c {output.gtf} > {output.gtf_gz}
        tabix -p gff -f {output.gtf_gz} ) &> {log}
        """

rule sqanti3:
    input:
        gtf1 = rules.collapse.output.gtf,
        gtf2 = config["HUMAN_ANNOTATION_GTF"],
        fa = config["HUMAN_GENOME_FASTA"]
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