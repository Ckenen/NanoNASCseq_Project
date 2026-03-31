#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
FQDIR = "results/2_demux/4_trimmed"
OUTDIR = "results/3_mapping"
# RUN_CELLS = RUN_CELLS[:1]
GROUPS1 = ["4_statClip", "5_markdup", "6_rmdup"]

rule all:
    input:
        expand(OUTDIR + "/0_index/{species}.mmi", species=["Human", "Mouse"]),
        expand(OUTDIR + "/{group1}/{run_cell}.bam", group1=GROUPS1, run_cell=RUN_CELLS),
        expand(OUTDIR + "/{group1}/{run_cell}.flagstat", group1=GROUPS1, run_cell=RUN_CELLS),

rule minimap2_build:
    input:
        fa = lambda wildcards: config["%s_GENOME_FASTA" % wildcards.species.upper()]
    output:
        mmi = OUTDIR + "/0_index/{species}.mmi"
    log:
        OUTDIR + "/0_index/{species}.log"
    threads:
        4
    shell:
        """
        minimap2 -t {threads} -x splice -d {output.mmi} {input.fa} &> {log}
        """

rule minimap2:
    input:
        fq = FQDIR + "/{run}/{cell}/trimmed.fastq.gz",
        mmi = lambda wildcards: OUTDIR + "/0_index/%s.mmi" % get_species(wildcards.cell),
        bed = lambda wildcards: config["%s_TRANSCRIPT_BED" % get_species(wildcards.cell).upper()]
    output:
        bam = OUTDIR + "/1_minimap2/{run}/{cell}.bam",
        bai = OUTDIR + "/1_minimap2/{run}/{cell}.bam.bai",
    log:
        OUTDIR + "/1_minimap2/{run}/{cell}.log"
    params:
        rg = '@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}'
    threads:
        12
    shell:
        """(
        minimap2 -ax splice -u f -Y --MD -R '{params.rg}' -t {threads} \
            --junc-bed {input.bed} {input.mmi} {input.fq} \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {output.bam} -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = OUTDIR + "/2_filtered/{run}/{cell}.bam",
        bai = OUTDIR + "/2_filtered/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/2_filtered/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            -q 30 -m 200 -F 2308 -o {output.bam} {input.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule extract_umi:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/3_extractUMI/{run}/{cell}.bam",
        bai = OUTDIR + "/3_extractUMI/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/3_extractUMI/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        ./scripts/mapping/extract_umi.py {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule stat_clip:
    input:
        bam = rules.extract_umi.output.bam
    output:
        bam = OUTDIR + "/4_statClip/{run}/{cell}.bam",
        bai = OUTDIR + "/4_statClip/{run}/{cell}.bam.bai",
        tsv = OUTDIR + "/4_statClip/{run}/{cell}.tsv"
    log:
        OUTDIR + "/4_statClip/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools StatClip -c 20 -s {output.tsv} -i {input.bam} -o {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_duplicate:  # Time-cost: 20221218_Blastocyst_69.C33
    input:
        bam = rules.stat_clip.output.bam
    output:
        bam = OUTDIR + "/5_markdup/{run}/{cell}.bam",
        bai = OUTDIR + "/5_markdup/{run}/{cell}.bam.bai",
        tsv = OUTDIR + "/5_markdup/{run}/{cell}.tsv"
    log:
        OUTDIR + "/5_markdup/{run}/{cell}.log"
    shell:
        """
        nasctools MarkDuplicate -s {output.tsv} -i {input.bam} -o {output.bam} &> {log}
        samtools index {output.bam}
        """

rule remove_duplicate:
    input:
        bam = rules.mark_duplicate.output.bam
    output:
        bam = OUTDIR + "/6_rmdup/{run}/{cell}.bam",
        bai = OUTDIR + "/6_rmdup/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/6_rmdup/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} -F 1024 -o {output.bam} {input.bam} &> {log}
        samtools index -@ {threads} {output.bam}
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
