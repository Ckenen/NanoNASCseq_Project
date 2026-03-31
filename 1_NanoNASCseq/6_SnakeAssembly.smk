#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
BAMDIR = "results/3_mapping/4_statClip"
OUTDIR = "results/6_assembly"
# RUN_CELLS = RUN_CELLS[:1]

rule all:
    input:
        expand(OUTDIR + "/1_stringtie/{run_cell}.gtf", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_sqanti3/{run_cell}", run_cell=RUN_CELLS),

rule stringtie:
    input:
        bam = BAMDIR + "/{run}/{cell}.bam",
        gtf = lambda wildcards: config["%s_ANNOTATION_GTF" % get_species(wildcards.cell).upper()]
    output:
        gtf = OUTDIR + "/1_stringtie/{run}/{cell}.gtf"
    log:
        OUTDIR + "/1_stringtie/{run}/{cell}.log"
    shell:
        """(
        stringtie {input.bam} -G {input.gtf} --fr -L \
            | sort -k1,1 -k4,4n | awk '$7!="."' > {output.gtf} ) &> {log}
        """

rule sqanti3:
    input:
        gtf1 = rules.stringtie.output.gtf,
        gtf2 = lambda wildcards: config["%s_ANNOTATION_GTF" % get_species(wildcards.cell).upper()],
        fasta = lambda wildcards: config["%s_GENOME_FASTA" % get_species(wildcards.cell).upper()]
    output:
        out = directory(OUTDIR + "/2_sqanti3/{run}/{cell}")
    log:
        OUTDIR + "/2_sqanti3/{run}/{cell}.log"
    conda:
        "SQANTI3.env"
    threads:
        4
    shell:
        """
        ../share/scripts/run_sqanti3_clean.sh {input} {threads} {output} &> {log}
        """
