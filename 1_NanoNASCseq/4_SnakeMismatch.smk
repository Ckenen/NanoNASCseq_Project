#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/3_mapping/5_mark_duplicate"
OUTDIR = "results/4_mismatch"
#RUN_CELLS = RUN_CELLS[:10]

rule all:
    input:
        expand(OUTDIR + "/1_tagged_events/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_ratio_all/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/3_ratio_rmdup/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/4_ratio_consensus/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/5_ratio_consensus_linkage/{run_cell}.tsv", run_cell=RUN_CELLS),

rule tag_events:
    input:
        bam = INDIR + "/{run}/{cell}.bam",
        bed = lambda wildcards: get_snp_bed(wildcards.cell)
    output:
        bam = OUTDIR + "/1_tagged_events/{run}/{cell}.bam",
        bai = OUTDIR + "/1_tagged_events/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/1_tagged_events/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools MarkEvent -t {threads} -s {input.bed} -i {input.bam} -o {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule report_mismatch_all:
    input:
        bam = rules.tag_events.output.bam
    output:
        tsv = OUTDIR + "/2_ratio_all/{run}/{cell}.tsv"
    log:
        OUTDIR + "/2_ratio_all/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools ReportMismatch -t {threads} -i {input.bam} -o {output.tsv} &> {log}
        """

rule report_mismatch_rmdup:
    input:
        bam = rules.tag_events.output.bam
    output:
        tsv = OUTDIR + "/3_ratio_rmdup/{run}/{cell}.tsv"
    log:
        OUTDIR + "/3_ratio_rmdup/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools ReportMismatch --discard-duplicates -t {threads} -i {input.bam} -o {output.tsv} &> {log}
        """

rule report_mismatch_consensus:
    input:
        bam = rules.tag_events.output.bam,
        fasta = lambda wildcards: get_genome_fasta(wildcards.cell)
    output:
        tsv = OUTDIR + "/4_ratio_consensus/{run}/{cell}.tsv",
        tsv2 = OUTDIR + "/4_ratio_consensus/{run}/{cell}.events.tsv"
    log:
        OUTDIR + "/4_ratio_consensus/{run}/{cell}.log"
    shell:
        """
        ./scripts/mismatch/report_mismatch_ratio.consensus.py {input.bam} {input.fasta} {output.tsv} {output.tsv2} &> {log}
        """

rule report_mismatch_consensus_linkage:
    input:
        tsv = rules.report_mismatch_consensus.output.tsv2
    output:
        tsv = OUTDIR + "/5_ratio_consensus_linkage/{run}/{cell}.tsv"
    log:
        OUTDIR + "/5_ratio_consensus_linkage/{run}/{cell}.log"
    shell:
        """
        ./scripts/mismatch/report_mismatch_ratio.consensus_linkage.py {input.tsv} {output.tsv} &> {log}
        """