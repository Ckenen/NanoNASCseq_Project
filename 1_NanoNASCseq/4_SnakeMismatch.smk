#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
BAMDIR = "results/3_mapping/5_markdup"
OUTDIR = "results/4_mismatch"
# RUN_CELLS = RUN_CELLS[:1]
GROUPS = ["1_all", "2_rmdup", "3_consensus", "4_consensus_linkage"]

rule all:
    input:
        expand(OUTDIR + "/1_tagged_events/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/1_tagged_events/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_mismatch_ratios/{group}/{run_cell}.tsv", group=GROUPS, run_cell=RUN_CELLS),

rule tag_events:
    input:
        bam = BAMDIR + "/{run}/{cell}.bam",
        bed = lambda wildcards: config["%s_SNP_BED_GZ" % get_species(wildcards.cell).upper()]
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
        tsv = OUTDIR + "/2_mismatch_ratios/1_all/{run}/{cell}.tsv"
    log:
        OUTDIR + "/2_mismatch_ratios/1_all/{run}/{cell}.log"
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
        tsv = OUTDIR + "/2_mismatch_ratios/2_rmdup/{run}/{cell}.tsv"
    log:
        OUTDIR + "/2_mismatch_ratios/2_rmdup/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools ReportMismatch --discard-duplicates -t {threads} -i {input.bam} -o {output.tsv} &> {log}
        """

rule report_mismatch_consensus:
    input:
        bam = rules.tag_events.output.bam,
        fasta = lambda wildcards: config["%s_GENOME_FASTA" % get_species(wildcards.cell).upper()]
    output:
        tsv = OUTDIR + "/2_mismatch_ratios/3_consensus/{run}/{cell}.tsv",
        tsv2 = OUTDIR + "/2_mismatch_ratios/3_consensus/{run}/{cell}.events.tsv"
    log:
        OUTDIR + "/2_mismatch_ratios/3_consensus/{run}/{cell}.log"
    shell:
        """
        ./scripts/mismatch/report_mismatch_ratio.consensus.py {input.bam} {input.fasta} 0.75 {output.tsv} {output.tsv2} &> {log}
        """

rule report_mismatch_consensus_linkage:
    input:
        tsv = rules.report_mismatch_consensus.output.tsv2
    output:
        tsv = OUTDIR + "/2_mismatch_ratios/4_consensus_linkage/{run}/{cell}.tsv"
    log:
        OUTDIR + "/2_mismatch_ratios/4_consensus_linkage/{run}/{cell}.log"
    shell:
        """
        ./scripts/mismatch/report_mismatch_ratio.consensus_linkage.py {input.tsv} {output.tsv} &> {log}
        """

rule flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """
