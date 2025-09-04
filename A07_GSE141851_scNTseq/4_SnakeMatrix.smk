#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
SAMPLES = config["SAMPLES"]
OUTDIR = "results/04_matrix"

rule all:
    input:
        expand(OUTDIR + "/01_TC/{sample}.tsv_q27.tsv_corrected.tsv", sample=SAMPLES),
        expand(OUTDIR + "/02_CorrectedBam/{sample}.TagTC.corrected.bam", sample=SAMPLES),
        expand(OUTDIR + "/03_Matrix/{sample}.TagTC.corrected_gene_cell_UMI_read.txt", sample=SAMPLES),
        expand(OUTDIR + "/03_Matrix/{sample}_TC_matrix.rds", sample=SAMPLES),
        expand(OUTDIR + "/04_h5ad/{sample}.h5ad", sample=SAMPLES),

def get_background_sample(sample):
    if "K562" in sample:
        return "SRR11683990_K562-4sU"
        # if "2ndSS" in sample:
        #     return "SRR11683993_K562-4sU-2ndSS"
        # else:
        #     return "SRR11683990_K562-4sU"
    elif "mESC" in sample:
        return "SRR12225122_mESC-4sU"
    elif "Mix" in sample:
        return "SRR10670112_Mix-4sU"
    elif "Cortical" in sample:
        return "SRR10670112_Mix-4sU"
    else:
        assert False

rule background_correction:
    input:
        bg = lambda wildcards: "results/03_mismatch/04_TC/%s.tsv_q27.tsv" % get_background_sample(wildcards.sample),
        tsv = "results/03_mismatch/04_TC/{sample}.tsv_q27.tsv"
    output:
        tsv1 = temp(OUTDIR + "/01_TC/{sample}.tsv_q27.tsv"),
        tsv2 = OUTDIR + "/01_TC/{sample}.tsv_q27.tsv_corrected.tsv"
    shell:
        """
        ln -s `readlink -f {input.tsv}` {output.tsv1}
        perl software/scNT-seq-master/TC_calling/scripts/background_correction.pl \
            -bg {input.bg} -in {output.tsv1}
        """

rule TagIntronicRead:
    input:
        tsv = rules.background_correction.output.tsv2,
        bam = "results/02_dropseq/12_DetectedErrors/{sample}.bam"
    output:
        bam = OUTDIR + "/02_CorrectedBam/{sample}.TagTC.corrected.bam"
    shell:
        """
        perl software/scNT-seq-master/TC_calling/scripts/TagIntronicRead_V5.pl \
            -read {input.tsv} -bam {input.bam}
        mv `basename {output.bam}` {output.bam}
        """

rule extract_digital_expression_matrix:
    input:
        bam = rules.TagIntronicRead.output.bam
    output:
        txt = OUTDIR + "/03_Matrix/{sample}.TagTC.corrected_gene_cell_UMI_read.txt"
    shell:
        """
        perl software/scNT-seq-master/TC_calling/scripts/extract_digital_expression_matrix.pl {input.bam}
        mv `basename {output.txt}` {output.txt}
        """

rule Generate_T_C_matrix:
    input:
        txt = rules.extract_digital_expression_matrix.output.txt
    output:
        rds = OUTDIR + "/03_Matrix/{sample}_TC_matrix.rds"
    shell:
        """
        Rscript software/scNT-seq-master/TC_calling/scripts/Generate_T_C_matrix.R {input.txt} 1200 {output.rds}
        """

rule make_h5ad:
    input:
        txt = rules.extract_digital_expression_matrix.output.txt
    output:
        h5ad = OUTDIR + "/04_h5ad/{sample}.h5ad"
    shell:
        """
        ./scripts/make_h5ad.py {input.txt} {output.h5ad}
        """