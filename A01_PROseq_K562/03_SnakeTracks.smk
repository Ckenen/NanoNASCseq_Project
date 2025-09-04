#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
OUTDIR = "results/03_tracks"

rule all:
    input:
        expand(OUTDIR + "/01_bigwigs_all/{sample}_raw_both.bw", sample=SAMPLES),
        expand(OUTDIR + "/02_bigwigs_rmdup/{sample}_raw_both.bw", sample=SAMPLES),

rule bam2bw_all:
    input:
        bam = "results/02_mapping/02_filtered/{sample}.bam"
    output:
        bw = OUTDIR + "/01_bigwigs_all/{sample}_raw_both.bw"
    params:
        prefix = OUTDIR + "/01_bigwigs_all/{sample}"
    shell:
        """
        ./scripts/bam2bw_proseq.sh {input.bam} {params.prefix}
        """

rule bam2bw_rmdup:
    input:
        bam = "results/02_mapping/03_rmdup/{sample}.bam"
    output:
        bw = OUTDIR + "/02_bigwigs_rmdup/{sample}_raw_both.bw"
    params:
        prefix = OUTDIR + "/02_bigwigs_rmdup/{sample}"
    shell:
        """
        ./scripts/bam2bw_proseq.sh {input.bam} {params.prefix}
        """
        