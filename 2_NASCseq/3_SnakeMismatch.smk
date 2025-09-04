#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
SNP_BED_GZ = config["SNP_BED_GZ"]
BAMDIR = "results/2_mapping/4_marked_strand"
OUTDIR = "results/3_mismatch"

rule all:
    input:
        expand(OUTDIR + "/1_marked_event/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_marked_new/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/3_mismatch_ratio/{run_cell}.tsv", run_cell=RUN_CELLS),

rule mark_events:
    input:
        bam = BAMDIR + "/{run}/{cell}.bam",
        bed = SNP_BED_GZ
    output:
        bam = OUTDIR + "/1_marked_event/{run}/{cell}.bam",
        bai = OUTDIR + "/1_marked_event/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/1_marked_event/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        nasctools MarkEvent -t {threads} -s {input.bed} -i {input.bam} -o {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_new:
    input:
        bam = rules.mark_events.output.bam,
    output:
        bam = OUTDIR + "/2_marked_new/{run}/{cell}.bam",
        bai = OUTDIR + "/2_marked_new/{run}/{cell}.bam.bai",
        tsv = OUTDIR + "/2_marked_new/{run}/{cell}.tsv"
    log:
        OUTDIR + "/2_marked_new/{run}/{cell}.log"
    shell:
        """
        nasctools MarkNew -s {output.tsv} -i {input.bam} -o {output.bam} &> {log}
        samtools index {output.bam}
        """

rule report_mismatch_ratio:
    input:
        bam = rules.mark_new.output.bam
    output:
        tsv = OUTDIR + "/3_mismatch_ratio/{run}/{cell}.tsv"
    log:
        OUTDIR + "/3_mismatch_ratio/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        nasctools ReportMismatch -t {threads} -i {input.bam} -o {output.tsv} &> {log}
        """
