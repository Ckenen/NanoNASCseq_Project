#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
BAMDIR = "results/3_mapping/2_filtered"
OUTDIR = "results/4_assembly"
# RUN_CELLS = RUN_CELLS[:1]

rule all:
    input:
        expand(OUTDIR + "/1_stringtie/{run_cell}.gtf.gz", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_sqanti3/{run_cell}", run_cell=RUN_CELLS),

rule stringtie:
    input:
        bam = BAMDIR + "/{run}/{cell}.bam",
        gtf = lambda wildcards: config["%s_ANNOTATION_GTF" % get_species(wildcards.cell).upper()]
    output:
        gtf1 = OUTDIR + "/1_stringtie/{run}/{cell}.gtf",
        gtf2 = OUTDIR + "/1_stringtie/{run}/{cell}.gtf.gz",
        tbi2 = OUTDIR + "/1_stringtie/{run}/{cell}.gtf.gz.tbi"
    log:
        OUTDIR + "/1_stringtie/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        stringtie {input.bam} -G {input.gtf} --fr -L | bedtools sort | awk '$7!="."' > {output.gtf1}
        bgzip -c {output.gtf1} > {output.gtf2}
        tabix -p gff {output.gtf2} ) &> {log}
        """

rule sqanti3:
    input:
        gtf_que = rules.stringtie.output.gtf1,
        gtf_ref = lambda wildcards: config["%s_ANNOTATION_GTF" % get_species(wildcards.cell).upper()],
        fasta = lambda wildcards: config["%s_GENOME_FASTA" % get_species(wildcards.cell).upper()]
    output:
        out = directory(OUTDIR + "/2_sqanti3/{run}/{cell}")
    log:
        OUTDIR + "/2_sqanti3/{run}/{cell}.log"
    threads:
        4
    conda:
        "SQANTI3.env"
    shell:
        """
        ../share/scripts/run_sqanti3_clean.sh {input.gtf_que} {input.gtf_ref} {input.fasta} {threads} {output.out} &> {log}
        """