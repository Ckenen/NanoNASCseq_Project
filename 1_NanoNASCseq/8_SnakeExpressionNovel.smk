#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"

EVENTDIR = "results/4_mismatch/2_mismatch_ratios/3_consensus"
ALLELEDIR = "results/5_expression/1_expressed_alleles"
BEDDIR = "results/5_expression/2_collapsed"

OUTDIR = "results/8_expression_novel"

# RUN_CELLS = RUN_CELLS[:1]
GROUPS = ['K562', 'mESC', 'MouseBlastocyst']
CUTOFFS = ["min_read_1_min_tc_1", "min_read_2_min_tc_1", "min_read_2_min_tc_2"]

RUN_CELLS = list(filter(lambda x: get_cell_type(x.split("/")[1]) in GROUPS, RUN_CELLS))

rule all:
    input:
        expand(OUTDIR + "/1_isoform_category/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_quant_isoforms/{cutoff}/{run_cell}.tsv", cutoff=CUTOFFS, run_cell=RUN_CELLS),

def get_novel_gtf(cell):
    ct = get_cell_type(cell)
    if ct == "K562" or ct == "K562_Mix" or ct == "K562_FUCCI":
        group = "K562"
    elif ct == "mESC" or ct == "mESC_Mix" or ct == "mESC_EXOSC2":
        group = "mESC"
    elif ct == "MouseBlastocyst" or ct == "Mouse2C" or ct == "MouseBigCell" or ct == "MouseMII":
        group = "MouseBlastocyst"
    else:
        print(ct)
        assert False
    return "results/7_assembly_novel/4_gtf/%s.all.gtf" % group

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
        event_tsv = EVENTDIR + "/{run}/{cell}.events.tsv",
        allele_tsv = ALLELEDIR +"/{run}/{cell}.tsv"
    output:
        tsv = OUTDIR + "/2_quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.tsv"
    log:
        OUTDIR + "/2_quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.log"
    shell:
        """
        ./scripts/expression/quant_isoforms.py {input.stat_tsv} {input.event_tsv} \
            {input.allele_tsv} {wildcards.size} {wildcards.tc} {output.tsv} &> {log}
        """
