#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
GTFDIR = "results/5_expression/2_collapsed"
BEDDIR = "results/5_expression/2_collapsed"
OUTDIR = "results/8_expression_novel"
RUN_CELLS = RUN_CELLS_CELLLINE + RUN_CELLS_BLASTOCYST

rule all:
    input:
        expand(OUTDIR + "/1_isoform_category/{run_cell}.tsv", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/2_quant_isoforms/min_read_1_min_tc_1/{run_cell}.tsv", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/2_quant_isoforms/min_read_2_min_tc_1/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_quant_isoforms/min_read_2_min_tc_2/{run_cell}.tsv", run_cell=RUN_CELLS),

def get_novel_gtf(cell):
    group = get_group(cell)
    return "results/7_assembly_custom/4_gtf/%s.all.gtf" % group

rule stat_isoform_category:
    input:
        gtf = lambda wildcards: get_novel_gtf(wildcards.cell),
        bed = BEDDIR + "/{run}/{cell}.bed.gz"
    output:
        tsv = OUTDIR + "/1_isoform_category/{run}/{cell}.tsv"
    log:
        OUTDIR + "/1_isoform_category/{run}/{cell}.log"
    shell:
        """
        ./scripts/expression/stat_isoform_category.py {input} {output} &> {log}
        """

rule quant_isoforms:
    input:
        stat_tsv = rules.stat_isoform_category.output.tsv,
        event_tsv = "results/4_mismatch/04_ratio_consensus/{run}/{cell}.events.tsv",
        allele_tsv = "results/5_expression/01_expressed_alleles/{run}/{cell}.tsv"
    output:
        tsv = OUTDIR + "/2_quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.tsv"
    log:
        OUTDIR + "/2_quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.log"
    shell:
        """
        ./scripts/expression/quant_isoforms.py {input.stat_tsv} {input.event_tsv} \
            {input.allele_tsv} {wildcards.size} {wildcards.tc} {output.tsv} &> {log}
        """
