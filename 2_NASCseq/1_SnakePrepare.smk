#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
BOWTIE2_INDEX = config["BOWTIE2_INDEX_RIBO"]
FQDIR = "data/datasets"
OUTDIR = "results/1_prepare"
# RUN_CELLS = RUN_CELLS[:4]

rule all:
    input:
        # expand(OUTDIR + "/1_cutadapt/{run_cell}_1.fastq.gz", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_bowtie2/{run_cell}.unmapped.fastq.1.gz", run_cell=RUN_CELLS),

rule cutadapt: # tanglab
    input:
        fq1 = FQDIR + "/{run}/{cell}_1.fastq.gz",
        fq2 = FQDIR + "/{run}/{cell}_2.fastq.gz"
    output:
        fq1 = OUTDIR + "/1_cutadapt/{run}/{cell}_1.fastq.gz",
        fq2 = OUTDIR + "/1_cutadapt/{run}/{cell}_2.fastq.gz"
    log:
        OUTDIR + "/1_cutadapt/{run}/{cell}.log"
    wildcard_constraints:
        run = "2022.+"
    conda:
        "cutadapt"
    threads:
        8
    shell:
        """
        cutadapt -j {threads} -q 30 -m 20 \
            -g GTGTATAAGAGACAG -g ATCAACGCAGAGTAC \
            -a CTGTCTCTTATACAC -a GTACTCTGCGTTGAT \
            -G GTGTATAAGAGACAG -G ATCAACGCAGAGTAC \
            -A CTGTCTCTTATACAC -A GTACTCTGCGTTGAT \
            -o {output.fq1} -p {output.fq2} \
            {input.fq1} {input.fq2} &> {log}
        """

rule cutadapt_GSE128273: # GSE128273_NASCseq and GSE128273_NASCseq_SE
    input:
        fq1 = FQDIR + "/{run}/{cell}_1.fastq.gz",
        fq2 = FQDIR + "/{run}/{cell}_2.fastq.gz"
    output:
        fq1 = temp(OUTDIR + "/1_cutadapt/{run}/{cell}_1.fastq.gz"),
        fq2 = temp(OUTDIR + "/1_cutadapt/{run}/{cell}_2.fastq.gz")
    log:
        OUTDIR + "/1_cutadapt/{run}/{cell}.log"
    wildcard_constraints:
        run = "GSE128273.+"
    conda:
        "cutadapt"
    threads:
        8
    shell:
        """
        cutadapt -j {threads} -q 30 -m 20 \
            -g GTGTATAAGAGACAG -g GCAGAGTACGGG \
            -a CTGTCTCTTATACAC -a CCCGTACTCTGC \
            -G GTGTATAAGAGACAG -G GCAGAGTACGGG \
            -A CTGTCTCTTATACAC -A CCCGTACTCTGC \
            -o {output.fq1} -p {output.fq2} \
            {input.fq1} {input.fq2} &> {log}
        """

rule bowtie2:
    input:
        fq1 = rules.cutadapt.output.fq1,
        fq2 = rules.cutadapt.output.fq2,
        idx = BOWTIE2_INDEX
    output:
        bam = OUTDIR + "/2_bowtie2/{run}/{cell}.bam",
        bai = OUTDIR + "/2_bowtie2/{run}/{cell}.bam.bai",
        fq1 = OUTDIR + "/2_bowtie2/{run}/{cell}.unmapped.fastq.1.gz",
        fq2 = OUTDIR + "/2_bowtie2/{run}/{cell}.unmapped.fastq.2.gz"
    log:
        OUTDIR + "/2_bowtie2/{run}/{cell}.log"
    conda:
        "bowtie2"
    params:
        prefix = OUTDIR + "/2_bowtie2/{run}/{cell}"
    threads:
        12
    shell:
        """
        bowtie2_keep_unmapped_fastq_pe.sh {input.fq1} {input.fq2} {input.idx} {threads} {params.prefix} &> {log}
        """
