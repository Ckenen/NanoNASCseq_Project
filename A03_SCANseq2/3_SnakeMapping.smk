#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
FQDIR = "results/2_demux/5_trim_polya"
OUTDIR = "results/3_mapping"
# RUN_CELLS = RUN_CELLS[:1]

rule all:
    input:
        # expand(OUTDIR + "/0_index/minimap2_{species}.mmi", species=["Human", "Mouse"]),
        # expand(OUTDIR + "/1_minimap2/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/1_minimap2/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_filtered/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_filtered/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/3_stat_clip/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/3_stat_clip/{run_cell}.flagstat", run_cell=RUN_CELLS),

rule minimap2_build:
    input:
        fa = lambda wildcards: config["%s_GENOME_FASTA" % wildcards.species.upper()]
    output:
        mmi = OUTDIR + "/0_index/minimap2_{species}.mmi"
    log:
        OUTDIR + "/0_index/minimap2_{species}.log"
    shell:
        """
        minimap2 -t 20 -x splice -d {output.mmi} {input.fa} &> {log}
        """    

rule minimap2:
    input:
        fq = FQDIR + "/{run}/{cell}.fastq.gz",
        mmi = lambda wildcards: OUTDIR + "/0_index/minimap2_%s.mmi" % get_species(wildcards.cell),
        bed = lambda wildcards: config["%s_TRANSCRIPT_BED" % get_species(wildcards.cell).upper()]
    output:
        bam = OUTDIR + "/1_minimap2/{run}/{cell}.bam",
        bai = OUTDIR + "/1_minimap2/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/1_minimap2/{run}/{cell}.log"
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
        bam = OUTDIR + "/2_filtered/{run}/{cell}.bam",
        bai = OUTDIR + "/2_filtered/{run}/{cell}.bam.bai"
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
        bam = OUTDIR + "/3_stat_clip/{run}/{cell}.bam",
        bai = OUTDIR + "/3_stat_clip/{run}/{cell}.bam.bai",
        txt = OUTDIR + "/3_stat_clip/{run}/{cell}.tsv"
    log:
        OUTDIR + "/3_stat_clip/{run}/{cell}.log"
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

        
