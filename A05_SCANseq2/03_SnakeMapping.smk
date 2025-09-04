#!/usr/bin/env runsnakemake
include: "00_SnakeCommon.smk"
FQDIR = "results/02_demux/05_trim_polya"
OUTDIR = "results/03_mapping"
# RUN_CELLS = RUN_CELLS[:1]

rule all:
    input:
        # expand(OUTDIR + "/01_minimap2/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/01_minimap2/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/02_filtered/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/02_filtered/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/03_stat_clip/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/03_stat_clip/{run_cell}.flagstat", run_cell=RUN_CELLS),

rule minimap2:
    input:
        fq = FQDIR + "/{run}/{cell}.fastq.gz",
        mmi = lambda wildcards: config["%s_SPLICE_MMI" % get_species(wildcards.cell).upper()],
        bed = lambda wildcards: config["%s_TRANSCRIPT_BED" % get_species(wildcards.cell).upper()]
    output:
        bam = temp(OUTDIR + "/01_minimap2/{run}/{cell}.bam"),
        bai = temp(OUTDIR + "/01_minimap2/{run}/{cell}.bam.bai")
    log:
        OUTDIR + "/01_minimap2/{run}/{cell}.log"
    params:
        rg = "@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}"
    threads:
        12
    shell:
        """
        minimap2 -a -x splice -u f -Y --MD -R "{params.rg}" \
            --junc-bed {input.bed} -t {threads} {input.mmi} {input.fq} \
            | samtools view --no-PG -@ {threads} -u - \
            | samtools sort --no-PG -@ {threads} -T {output.bam} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = OUTDIR + "/02_filtered/{run}/{cell}.bam",
        bai = OUTDIR + "/02_filtered/{run}/{cell}.bam.bai"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            -q 30 -m 200 -F 2308 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule stat_clip:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/03_stat_clip/{run}/{cell}.bam",
        bai = OUTDIR + "/03_stat_clip/{run}/{cell}.bam.bai",
        txt = OUTDIR + "/03_stat_clip/{run}/{cell}.tsv"
    log:
        OUTDIR + "/03_stat_clip/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools StatClip -c 20 -s {output.txt} -i {input.bam} -o {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule bam_flagstat:
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

        
