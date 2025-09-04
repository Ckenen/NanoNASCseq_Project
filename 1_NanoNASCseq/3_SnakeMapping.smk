#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
FQDIR = "results/2_demux/4_trimmed"
OUTDIR = "results/3_mapping"

rule all:
    input:
        # expand(OUTDIR + "/1_minimap2/{run_cell}.bam", run_cell=RUN_CELLS),
         expand(OUTDIR + "/1_minimap2/{run_cell}.flagstat", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/2_filtered/{run_cell}.bam", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/3_extract_umi/{run_cell}.bam", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/4_stat_clip/{run_cell}.bam", run_cell=RUN_CELLS),
        #expand(OUTDIR + "/5_mark_duplicate/{run_cell}.bam", run_cell=RUN_CELLS),
        #expand(OUTDIR + "/5_mark_duplicate/{run_cell}.flagstat", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/6_remove_duplicate/{run_cell}.bam", run_cell=RUN_CELLS),

rule minimap2:
    input:
        fq = FQDIR + "/{run}/{cell}/trimmed.fastq.gz",
        mmi = lambda wildcards: get_genome_splice_mmi(wildcards.cell),
        bed = lambda wildcards: get_transcript_bed(wildcards.cell)
    output:
        bam = temp(OUTDIR + "/1_minimap2/{run}/{cell}.bam"),
        bai = temp(OUTDIR + "/1_minimap2/{run}/{cell}.bam.bai"),
        txt = OUTDIR + "/1_minimap2/{run}/{cell}.flagstat"
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
        samtools index -@ {threads} {output.bam} 
        samtools flagstat -@ {threads} {output.bam} > {output.txt}  ) &> {log}
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
        samtools view -@ {threads} --expr 'rname =~ "^(chr([0-9]+|[XY])|FUCCI|ctgEXOSC2)$"' \
            -q 30 -m 200 -F 2308 -o {output.bam} {input.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule extract_umi:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/3_extract_umi/{run}/{cell}.bam",
        bai = OUTDIR + "/3_extract_umi/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/3_extract_umi/{run}/{cell}.log"
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
        bam = OUTDIR + "/4_stat_clip/{run}/{cell}.bam",
        bai = OUTDIR + "/4_stat_clip/{run}/{cell}.bam.bai",
        tsv = OUTDIR + "/4_stat_clip/{run}/{cell}.tsv"
    log:
        OUTDIR + "/4_stat_clip/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools StatClip -c 20 -s {output.tsv} -i {input.bam} -o {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

# Time-cost: 20221218_Blastocyst_69.C33

rule mark_duplicate: 
    input:
        bam = rules.stat_clip.output.bam
    output:
        bam = OUTDIR + "/5_mark_duplicate/{run}/{cell}.bam",
        bai = OUTDIR + "/5_mark_duplicate/{run}/{cell}.bam.bai",
        tsv = OUTDIR + "/5_mark_duplicate/{run}/{cell}.tsv"
    log:
        OUTDIR + "/5_mark_duplicate/{run}/{cell}.log"
    shell:
        """
        nasctools MarkDuplicate -s {output.tsv} -i {input.bam} -o {output.bam} &> {log}
        samtools index {output.bam}
        """

rule remove_duplicate:
    input:
        bam = rules.mark_duplicate.output.bam
    output:
        bam = OUTDIR + "/6_remove_duplicate/{run}/{cell}.bam",
        bai = OUTDIR + "/6_remove_duplicate/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/6_remove_duplicate/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} -F 1024 -o {output.bam} {input.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

# rule flagstat:
#     input:
#         bam = "{prefix}.bam"
#     output:
#         txt = "{prefix}.flagstat"
#     threads:
#         4
#     shell:
#         """
#         samtools flagstat -@ {threads} {input.bam} > {output.txt}
#         """
