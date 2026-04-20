#!/usr/bin/env runsnakemake
configfile: "config.yaml"
import pandas as pd
RUNS = config["RUNS"]
DAT = pd.read_excel("data/NASCseq.xlsx")
DAT = DAT[DAT["Run"].isin(RUNS)]
DAT["RunCell"] = ["%s/%s" % (run, cell) for run, cell in DAT[["Run", "Cell"]].values]
RUN_CELLS = DAT["RunCell"]
GROUPS = ["filtered", "rmdup", "markstrand", "markevent", "marknew"]
OUTDIR = "results/1_nascseq"
# RUN_CELLS = RUN_CELLS[:1]

rule all:
    input:
        OUTDIR + "/0_index/ribosomal_bowtie2.idx",
        OUTDIR + "/0_index/genome_star.idx",
        # expand(OUTDIR + "/1_fastqs/cutadapt/{run_cell}_1.fastq.gz", run_cell=RUN_CELLS),
        expand(OUTDIR + "/1_fastqs/bowtie2/{run_cell}.unmapped.fastq.1.gz", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_bams/star/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_bams/{A}/{run_cell}.bam", A=GROUPS, run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_bams/{A}/{run_cell}.flagstat", run_cell=RUN_CELLS, A=GROUPS),
        expand(OUTDIR + "/3_mismatch_ratio/marknew/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/4_quantify/fpkm/marknew/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/5_stats/gene_number/marknew/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/5_stats/read_location/marknew/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/5_stats/snr/marknew/{run_cell}.csv", run_cell=RUN_CELLS),

# Build index

rule bowtie2_build:
    input:
        fa = config["FASTA_RIBO"]
    output:
        idx = directory(OUTDIR + "/0_index/ribosomal_bowtie2.idx")
    log:
        OUTDIR + "/0_index/ribosomal_bowtie2.idx.log"
    conda:
        "bowtie2"
    shell:
        """(
        mkdir {output.idx}
        bowtie2-build {input.fa} {output.idx}/ref ) &> {log}
        """

rule star_build:
    input:
        fa = config["FASTA"],
        gtf = config["GTF"]
    output:
        idx = directory(OUTDIR + "/0_index/genome_star.idx")
    log:
        OUTDIR + "/0_index/genome_star.log"
    conda:
        "star"
    threads:
        24
    shell:
        """
        mkdir -p {output.idx}
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.idx} \
            --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} &> {log}
        """

# Pipeline

rule cutadapt: # tanglab datasets
    input:
        fq1 = config["FQDIR"] + "/{run}/{cell}_1.fastq.gz",
        fq2 = config["FQDIR"] + "/{run}/{cell}_2.fastq.gz"
    output:
        fq1 = OUTDIR + "/1_fastqs/cutadapt/{run}/{cell}_1.fastq.gz",
        fq2 = OUTDIR + "/1_fastqs/cutadapt/{run}/{cell}_2.fastq.gz"
    log:
        OUTDIR + "/1_fastqs/cutadapt/{run}/{cell}.log"
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
        fq1 = config["FQDIR"] + "/{run}/{cell}_1.fastq.gz",
        fq2 = config["FQDIR"] + "/{run}/{cell}_2.fastq.gz"
    output:
        fq1 = OUTDIR + "/1_fastqs/cutadapt/{run}/{cell}_1.fastq.gz",
        fq2 = OUTDIR + "/1_fastqs/cutadapt/{run}/{cell}_2.fastq.gz"
    log:
        OUTDIR + "/1_fastqs/cutadapt/{run}/{cell}.log"
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
        idx = rules.bowtie2_build.output.idx
    output:
        bam = temp(OUTDIR + "/1_fastqs/bowtie2/{run}/{cell}.bam"),
        bai = temp(OUTDIR + "/1_fastqs/bowtie2/{run}/{cell}.bam.bai"),
        fq1 = OUTDIR + "/1_fastqs/bowtie2/{run}/{cell}.unmapped.fastq.1.gz",
        fq2 = OUTDIR + "/1_fastqs/bowtie2/{run}/{cell}.unmapped.fastq.2.gz"
    log:
        OUTDIR + "/1_fastqs/bowtie2/{run}/{cell}.log"
    conda:
        "bowtie2"
    params:
        prefix = OUTDIR + "/1_fastqs/bowtie2/{run}/{cell}"
    threads:
        12
    shell:
        """
        ../share/scripts/bowtie2_keep_unmapped_fastq_pe.sh {input.fq1} {input.fq2} {input.idx} {threads} {params.prefix} &> {log}
        """

rule star:
    input:
        fq1 = rules.bowtie2.output.fq1,
        fq2 = rules.bowtie2.output.fq2,
        idx = rules.star_build.output.idx
    output:
        bam = OUTDIR + "/2_bams/star/{run}/{cell}.bam"
    log:
        OUTDIR + "/2_bams/star/{run}/{cell}.log"
    params:
        prefix = OUTDIR + "/2_bams/star/{run}/{cell}"
    conda:
        "star"
    threads:
        12
    shell:
        """(
        STAR --runThreadN {threads} \
            --outFileNamePrefix {params.prefix}. \
            --genomeDir {input.idx} \
            --readFilesCommand zcat \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fq1} {input.fq2}
        mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.star.output.bam
    output:
        bam = OUTDIR + "/2_bams/filtered/{run}/{cell}.bam",
        bai = OUTDIR + "/2_bams/filtered/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/2_bams/filtered/{run}/{cell}.log"
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
        bam = OUTDIR + "/2_bams/rmdup/{run}/{cell}.bam",
        bai = OUTDIR + "/2_bams/rmdup/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/2_bams/rmdup/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sambamba markdup -t {threads} -r {input.bam} {output.bam} &> {log}
        """

rule mark_strand:
    input:
        bam = rules.rmdup.output.bam,
        bed = config["TRANSCRIPT_BED_GZ"]
    output:
        bam = OUTDIR + "/2_bams/markstrand/{run}/{cell}.bam",
        bai = OUTDIR + "/2_bams/markstrand/{run}/{cell}.bam.bai",
        tsv = OUTDIR + "/2_bams/markstrand/{run}/{cell}.tsv"
    log:
        OUTDIR + "/2_bams/markstrand/{run}/{cell}.log"
    shell:
        """
        nasctools MarkStrand -s U -m {output.tsv} -g {input.bed} -i {input.bam} -o {output.bam} &> {log}
        samtools index {output.bam}
        """

rule mark_events:
    input:
        bam = rules.mark_strand.output.bam,
        bed = config["SNP_BED_GZ"]
    output:
        bam = OUTDIR + "/2_bams/markevent/{run}/{cell}.bam",
        bai = OUTDIR + "/2_bams/markevent/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/2_bams/markevent/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        nasctools MarkEvent -t {threads} -s {input.bed} -i {input.bam} -o {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_new:
    input:
        bam = rules.mark_events.output.bam,
    output:
        bam = OUTDIR + "/2_bams/marknew/{run}/{cell}.bam",
        bai = OUTDIR + "/2_bams/marknew/{run}/{cell}.bam.bai",
        tsv = OUTDIR + "/2_bams/marknew/{run}/{cell}.tsv"
    log:
        OUTDIR + "/2_bams/marknew/{run}/{cell}.log"
    shell:
        """
        nasctools MarkNew -s {output.tsv} -i {input.bam} -o {output.bam} &> {log}
        samtools index {output.bam}
        """

# Common rules

rule report_mismatch_ratio:
    input:
        bam = OUTDIR + "/2_bams/{A}/{run}/{cell}.bam",
    output:
        tsv = OUTDIR + "/3_mismatch_ratio/{A}/{run}/{cell}.tsv"
    log:
        OUTDIR + "/3_mismatch_ratio/{A}/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        nasctools ReportMismatch -t {threads} -i {input.bam} -o {output.tsv} &> {log}
        """

rule calculate_fpkm:
    input:
        bam = OUTDIR + "/2_bams/{A}/{run}/{cell}.bam",
        bed = config["TRANSCRIPT_BED_GZ"],
        tsv = config["TRANSCRIPT_INFO"]
    output:
        txt = OUTDIR + "/4_quantify/fpkm/{A}/{run}/{cell}.tsv",
    log:
        OUTDIR + "/4_quantify/fpkm/{A}/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        calculate_fpkm.py --threads {threads} --annotation {input.tsv} \
            --stranded TAG --strand-tag ST --new \
            {input.bam} {input.bed} {output.txt} &> {log}
        """

rule report_gene_number:
    input:
        tsv = OUTDIR + "/4_quantify/fpkm/{A}/{run}/{cell}.tsv",
    output:
        tsv = OUTDIR + "/5_stats/gene_number/{A}/{run}/{cell}.tsv"
    shell:
        """
        ../share/scripts/report_gene_number.py {input} > {output}
        """

rule stat_read_location:
    input:
        bam = OUTDIR + "/2_bams/{A}/{run}/{cell}.bam",
        bed = config["TRANSCRIPT_BED_GZ"]
    output:
        bed1 = temp(OUTDIR + "/5_stats/read_location/{A}/{run}/{cell}.bed"),
        bed2 = OUTDIR + "/5_stats/read_location/{A}/{run}/{cell}.bed.gz",
        tsv = OUTDIR + "/5_stats/read_location/{A}/{run}/{cell}.tsv"
    log:
        OUTDIR + "/5_stats/read_location/{A}/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        stat_read_location.py -b {input.bam} -g {input.bed} -t {threads} -o {output.bed1} -s {output.tsv} &> {log}
        bgzip -c {output.bed1} > {output.bed2}
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

def get_pe_model(cell):
    if cell.startswith("GSE128273"):
        return "../0_Analysis/4_signal_to_noise/results/models/GSE128273_NASCseq.K562.s4U_0uM_3h.pkl"
    else:
        return "../0_Analysis/4_signal_to_noise/results/models/NASCseq.K562.s4U_0uM_3h.pkl"

rule stat_snr:
    input:
        bam = OUTDIR + "/2_bams/{A}/{run}/{cell}.bam",
        ratio = OUTDIR + "/3_mismatch_ratio/{A}/{run}/{cell}.tsv",
        model = lambda wildcards: get_pe_model(wildcards.cell)
    output:
        ratio = OUTDIR + "/5_stats/snr/{A}/{run}/{cell}.ratio.csv",
        event = OUTDIR + "/5_stats/snr/{A}/{run}/{cell}.event.csv",
        tsv = OUTDIR + "/5_stats/snr/{A}/{run}/{cell}.csv"
    log:
        OUTDIR + "/5_stats/snr/{A}/{run}/{cell}.log"
    shell:
        """(
        ../share/scripts/make_one_row_mismatch_ratio.py -i {input.ratio} -o {output.ratio}
        ../share/scripts/make_t_tc_count.ngs.py {input.bam} {output.event}
        nasctools EstimateSNR -r {output.ratio} -e {output.event} -m {input.model} -o {output.tsv} ) &> {log}
        """
