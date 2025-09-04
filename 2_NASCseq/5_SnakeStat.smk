#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
OUTDIR = "results/5_stat"
# RUN_CELLS = RUN_CELLS[:1]

rule all:
    input:
        expand(OUTDIR + "/1_detected_genes/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_read_location/{run_cell}_summary.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/3_snr/{run_cell}.snr.csv", run_cell=RUN_CELLS),

rule report_gene_number:
    input:
        tsv = "results/04_expression/01_fpkm/{run}/{cell}.tsv"
    output:
        tsv = OUTDIR + "/1_detected_genes/{run}/{cell}.tsv"
    shell:
        """
        ./scripts/report_gene_number.py {input} > {output}
        """

rule stat_read_location:
    input:
        bam = "esults/3_mismatch/2_marked_new/{run}/{cell}.bam",
        bed = config["TRANSCRIPT_BED_GZ"]
    output:
        bed1 = temp(OUTDIR + "/2_read_location/{run}/{cell}.bed"),
        bed2 = OUTDIR + "/2_read_location/{run}/{cell}.bed.gz",
        tsv = OUTDIR + "/2_read_location/{run}/{cell}_summary.tsv",
    log:
        OUTDIR + "/2_read_location/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        stat_read_location.py -b {input.bam} -g {input.bed} -t {threads} -o {output.bed1} -s {output.tsv} &> {log}
        bgzip -c {output.bed1} > {output.bed2}
        """

def get_pe_model(cell):
    if cell.startswith("GSE128273"):
        return "../0_Analysis/4_signal_to_noise/results/models/GSE128273_NASCseq.K562.s4U_0uM_3h.pkl"
    else:
        return "../0_Analysis/4_signal_to_noise/results/models/NASCseq.K562.s4U_0uM_3h.pkl"

rule stat_snr:
    input:
        ratio = "results/3_mismatch/3_mismatch_ratio/{run}/{cell}.tsv",
        bam = "results/3_mismatch/2_marked_new/{run}/{cell}.bam",
        model = lambda wildcards: get_pe_model(wildcards.cell)
    output:
        ratio = OUTDIR + "/3_snr/{run}/{cell}.ratio.csv",
        event = OUTDIR + "/3_snr/{run}/{cell}.event.csv",
        tsv = OUTDIR + "/3_snr/{run}/{cell}.snr.csv"
    log:
        OUTDIR + "/3_snr/{run}/{cell}.log"
    shell:
        """(
        ../1_NanoNASCseq/scripts/stat/make_one_row_mismatch_ratio.py -i {input.ratio} -o {output.ratio}
        ../1_NanoNASCseq/scripts/stat/make_t_tc_count_ngs.py {input.bam} {output.event}
        nasctools EstimateSNR -r {output.ratio} -e {output.event} -m {input.model} -o {output.tsv} ) &> {log}
        """