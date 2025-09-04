#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
OUTDIR = "results/9_stat"
# RUN_CELLS = RUN_CELLS_K562 + RUN_CELLS_MESC
# RUN_CELLS = RUN_CELLS[:20]

rule all:
    input:
        #expand(OUTDIR + "/1_read_location/{run_cell}_summary.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_snr_raw/{run_cell}.snr.csv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_snr_corrected/{run_cell}.snr.csv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_snr_linkage/{run_cell}.snr.csv", run_cell=RUN_CELLS),
        #expand(OUTDIR + "/3_depth/{run_cell}.tsv", run_cell=RUN_CELLS),
        #expand(OUTDIR + "/4_chrom_reads/{run_cell}.tsv", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/5_consensus_accuracy/{run_cell}", run_cell=list(filter(lambda x: x.startswith("20220719_K562_3"), RUN_CELLS))),


rule stat_read_location:
    input:
        rd = "results/5_expression/2_collapsed/{run}/{cell}.bed.gz",
        tsv = "results/4_mismatch/4_ratio_consensus/{run}/{cell}.events.tsv",
        bed = lambda wildcards: config["%s_TRANSCRIPT_BED_GZ" % get_species(wildcards.cell).upper()]
    output:
        bed1 = temp(OUTDIR + "/1_read_location/{run}/{cell}.bed"),
        bed2 = OUTDIR + "/1_read_location/{run}/{cell}.bed.gz",
        tsv = OUTDIR + "/1_read_location/{run}/{cell}_summary.tsv",
    log:
        OUTDIR + "/1_read_location/{run}/{cell}.log"
    shell:
        """
        stat_read_location.py -l {input.rd} -e {input.tsv} -g {input.bed} -o {output.bed1} -s {output.tsv} &> {log}
        bgzip -c {output.bed1} > {output.bed2}
        """

def get_pe_model(cell, level):
    if get_species(cell) == "Human":
        return "../0_Analysis/4_signal_to_noise/results/models/NanoNASCseq.K562.s4U_0uM_3h.%s.pkl" % level
    else:
        return "../0_Analysis/4_signal_to_noise/results/models/NanoNASCseq.mESC.s4U_0uM_3h.%s.pkl" % level

rule stat_snr_raw:
    input:
        ratio = "results/4_mismatch/3_ratio_rmdup/{run}/{cell}.tsv",
        bam = "results/4_mismatch/1_tagged_events/{run}/{cell}.bam",
        model = lambda wildcards: get_pe_model(wildcards.cell, "raw"),
    output:
        ratio = OUTDIR + "/2_snr_raw/{run}/{cell}.ratio.csv",
        event = OUTDIR + "/2_snr_raw/{run}/{cell}.event.csv",
        tsv = OUTDIR + "/2_snr_raw/{run}/{cell}.snr.csv"
    log:
        OUTDIR + "/2_snr_raw/{run}/{cell}.log"
    shell:
        """(
        ./scripts/stat/make_one_row_mismatch_ratio.py -i {input.ratio} -o {output.ratio}
        ./scripts/stat/make_t_tc_count_ngs.py {input.bam} {output.event}
        nasctools EstimateSNR -r {output.ratio} -e {output.event} -m {input.model} -o {output.tsv} ) &> {log}
        """

rule stat_snr_corrected:
    input:
        ratio = "results/4_mismatch/4_ratio_consensus/{run}/{cell}.tsv",
        event = "results/4_mismatch/4_ratio_consensus/{run}/{cell}.events.tsv",
        model = lambda wildcards: get_pe_model(wildcards.cell, "corrected"),
    output:
        ratio = OUTDIR + "/2_snr_corrected/{run}/{cell}.ratio.csv",
        event = OUTDIR + "/2_snr_corrected/{run}/{cell}.event.csv",
        tsv = OUTDIR + "/2_snr_corrected/{run}/{cell}.snr.csv",
    log:
        OUTDIR + "/2_snr_corrected/{run}/{cell}.log"
    shell:
        """(
        ./scripts/stat/make_one_row_mismatch_ratio.py -i {input.ratio} -o {output.ratio}
        ./scripts/stat/make_t_tc_count_tgs.py {input.event} {output.event}
        nasctools EstimateSNR -r {output.ratio} -e {output.event} -m {input.model} -o {output.tsv} ) &> {log}
        """

rule stat_snr_linkage:
    input:
        ratio = "results/4_mismatch/5_ratio_consensus_linkage/{run}/{cell}.tsv",
        event = "results/4_mismatch/4_ratio_consensus/{run}/{cell}.events.tsv",
        model = lambda wildcards: get_pe_model(wildcards.cell, "linkage"),
    output:
        ratio = OUTDIR + "/2_snr_linkage/{run}/{cell}.ratio.csv",
        event = OUTDIR + "/2_snr_linkage/{run}/{cell}.event.csv",
        tsv = OUTDIR + "/2_snr_linkage/{run}/{cell}.snr.csv",
    log:
        OUTDIR + "/2_snr_linkage/{run}/{cell}.log"
    shell:
        """(
        ./scripts/stat/make_one_row_mismatch_ratio.py -i {input.ratio} -o {output.ratio}
        ./scripts/stat/make_t_tc_count_tgs.linkage.py {input.event} {output.event}
        nasctools EstimateSNR -r {output.ratio} -e {output.event} -m {input.model} -o {output.tsv} ) &> {log}
        """

rule stat_fq_depth:
    input:
        fq = "results/2_demux/4_trimmed/{run}/{cell}/trimmed.fastq.gz"
    output:
        tsv = OUTDIR + "/3_depth/{run}/{cell}.tsv"
    shell:
        """
        ./scripts/stat/stat_fq_depth.sh {input.fq} > {output.tsv}
        """

rule stat_chrom_reads:
    input:
        bam = "results/3_mapping/1_minimap2/{run}/{cell}.bam"
    output:
        tsv = OUTDIR + "/4_chrom_reads/{run}/{cell}.tsv"
    threads:
        4
    shell:
        """
        ./scripts/stat/stat_chrom_reads.sh {input.bam} {output.tsv}
        """

rule estimate_accuracy:
    input:
        bam = "results/3_mapping/5_mark_duplicate/{run}/{cell}.bam",
        mmi = lambda wildcards: get_genome_splice_mmi(wildcards.cell),
        bed = lambda wildcards: get_transcript_bed(wildcards.cell)
    output:
        out = directory(OUTDIR + "/5_consensus_accuracy/{run}/{cell}")
    log:
        OUTDIR + "/5_consensus_accuracy/{run}/{cell}.log"
    conda:
        "consensus"
    threads:
        24
    shell:
        """
        ./scripts/stat/estimate_accuracy.py \
            -b {input.bam} -m {input.mmi} \
            -g {input.bed} -t {threads} \
            -o {output} &> {log}
        """
