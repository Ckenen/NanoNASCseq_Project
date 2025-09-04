#!/usr/bin/env runsnakemake
SAMPLES = [
    "SRR11683990_K562-4sU",
    "SRR11683991_K562-TFEA",
    "SRR11683992_K562-4sU-TFEA",
    # "SRR11683993_K562-4sU-2ndSS",
    # "SRR11683994_K562-TFEA-2ndSS",
    # "SRR11683995_K562-4sU-TFEA-2ndSS",
]
OUTDIR = "results/signal_to_noise"
SAMPLE_CELLS = []
for s in SAMPLES:
    path = OUTDIR + "/split/%s/barcodes.txt" % s
    # print(path)
    if os.path.exists(path):
        for line in open(path):
            SAMPLE_CELLS.append("%s/%s" % (s, line.strip("\n").split("\t")[0]))
print("Cells:", len(SAMPLE_CELLS))

rule all:
    input:
        expand(OUTDIR + "/events/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/split/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/ratio/{sample_cell}.tsv", sample_cell=SAMPLE_CELLS),
        expand(OUTDIR + "/pc/{sample_cell}.tsv", sample_cell=SAMPLE_CELLS),

rule get_events:
    input:
        bam = "results/dropseq/12_DetectedErrors/{sample}.bam",
        bed = "/home/chenzonggui/species/homo_sapiens/hg38/snp151.3.lite.bed.gz"
    output:
        bam = OUTDIR + "/events/{sample}.bam"
    log:
        OUTDIR + "/events/{sample}.log"
    threads:
        12
    shell:
        """
        nasctools GetEvent -t {threads} -s {input.bed} {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule split_bam:
    input:
        bam = rules.get_events.output.bam
    output:
        out = directory(OUTDIR + "/split/{sample}")
    log:
        OUTDIR + "/split/{sample}.log"
    shell:
        """
        ./scripts/split_single_cell_bam.py {input.bam} {output} &> {log}
        """

rule report_ratio:
    input:
        bamdir = rules.split_bam.output.out
    output:
        tsv = OUTDIR + "/ratio/{sample}/{cell}.tsv"
    log:
        OUTDIR + "/ratio/{sample}/{cell}.log"
    threads:
        1
    shell:
        """
        nasctools ReportMismatch -t {threads} -s F {input.bamdir}/{wildcards.cell}.bam {output.tsv} &> {log}
        """

rule estimate_pc:
    input:
        tsv1 = "reports/Estimate.Pe.Model.scNTseq.tsv",
        tsv2 = rules.report_ratio.output.tsv,
        bamdir = rules.split_bam.output.out
    output:
        txt = OUTDIR + "/pc/{sample}/{cell}.tsv"
    shell:
        """
        ./scripts/estimate_pc.py -m short -b {input.bamdir}/{wildcards.cell}.bam {input.tsv1} {input.tsv2} > {output}
        """
