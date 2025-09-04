#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
FQDIR = "results/1_prepare/2_bowtie2"
STAR_GENOME = config["STAR_GENOME"]
TRANSCRIPT_BED_GZ = config["TRANSCRIPT_BED_GZ"]
OUTDIR = "results/2_mapping"
# RUN_CELLS = RUN_CELLS[:10]

rule all:
    input:
        # expand(OUTDIR + "/1_star_mapped/{run_cell}", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/2_bam_filtered/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_bam_filtered/{run_cell}.flagstat", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/3_rmdup/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/4_marked_strand/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/4_marked_strand/{run_cell}.flagstat", run_cell=RUN_CELLS),

rule star_mapping:
    input:
        fq1 = FQDIR + "/{run}/{cell}.unmapped.fastq.1.gz",
        fq2 = FQDIR + "/{run}/{cell}.unmapped.fastq.2.gz",
        idx = STAR_GENOME
    output:
        directory(OUTDIR + "/1_star_mapped/{run}/{cell}")
    conda:
        "star"
    log:
        OUTDIR + "/1_star_mapped/{run}/{cell}.log"
    threads:
        12
    shell:
        """
        mkdir {output}
        STAR --runThreadN {threads} \
            --outFileNamePrefix {output}/{wildcards.cell}. \
            --genomeDir {input.idx} \
            --genomeLoad LoadAndKeep \
            --readFilesCommand zcat \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

rule filter_bam:
    input:
        rules.star_mapping.output
    output:
        bam = OUTDIR + "/2_bam_filtered/{run}/{cell}.bam",
        bai = OUTDIR + "/2_bam_filtered/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/2_bam_filtered/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} -q 30 -d "NH:1" \
            --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            -f 2 -F 2308 -o {output.bam} \
            {input}/{wildcards.cell}.Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

# This step will affect final results. Because the uniq read are 
# randomly selected and do not consider the mismatch events.

rule rmdup: 
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/3_rmdup/{run}/{cell}.bam",
        bai = OUTDIR + "/3_rmdup/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/3_rmdup/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sambamba markdup -t {threads} -r {input.bam} {output.bam} &> {log}
        """

# rule mark_strand:
#     input:
#         bam = rules.rmdup.output.bam,
#         bed = TRANSCRIPT_BED_GZ
#     output:
#         bam = OUTDIR + "/4_marked_strand/{run}/{cell}.bam",
#         bai = OUTDIR + "/4_marked_strand/{run}/{cell}.bam.bai",
#         tsv = OUTDIR + "/4_marked_strand/{run}/{cell}.tsv"
#     log:
#         OUTDIR + "/4_marked_strand/{run}/{cell}.log"
#     shell:
#         """
#         nasctools MarkStrand -s U -r -m {output.tsv} -g {input.bed} -i {input.bam} -o {output.bam} &> {log}
#         samtools index {output.bam}
#         """

rule mark_strand:
    input:
        bam = rules.rmdup.output.bam,
        bed = TRANSCRIPT_BED_GZ
    output:
        bam = OUTDIR + "/4_marked_strand/{run}/{cell}.bam",
        bai = OUTDIR + "/4_marked_strand/{run}/{cell}.bam.bai",
        tsv = OUTDIR + "/4_marked_strand/{run}/{cell}.tsv"
    log:
        OUTDIR + "/4_marked_strand/{run}/{cell}.log"
    shell:
        """
        nasctools MarkStrand -s U -m {output.tsv} -g {input.bed} -i {input.bam} -o {output.bam} &> {log}
        samtools index {output.bam}
        """

# Common rules

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
