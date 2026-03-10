#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
BAMDIR = "results/3_mapping/3_bam_filtered"
OUTDIR = "results/5_stat"

rule all:
    input:
        expand(OUTDIR + "/1_read_location/{run_cell}_summary.tsv", run_cell=RUN_CELLS),

rule stat_read_location:
    input:
        bam = BAMDIR + "/{run}/{cell}.bam",
        bed = lambda wildcards: config["%s_TRANSCRIPT_BED_GZ" % get_species(wildcards.cell).upper()]
    output:
        bed1 = temp(OUTDIR + "/1_read_location/{run}/{cell}.bed"),
        bed2 = OUTDIR + "/1_read_location/{run}/{cell}.bed.gz",
        tbi2 = OUTDIR + "/1_read_location/{run}/{cell}.bed.gz.tbi",
        tsv = OUTDIR + "/1_read_location/{run}/{cell}_summary.tsv",
    log:
        OUTDIR + "/1_read_location/{run}/{cell}.log"
    shell:
        """
        stat_read_location.py -l {input.bam} -g {input.bed} -o {output.bed1} -s {output.tsv} &> {log}
        bgzip -c {output.bed1} > {output.bed2}
        """