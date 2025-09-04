#!/usr/bin/env runsnakemake
include: "00_SnakeCommon.smk"
RUNS = config["SAMPLES"]
FQDIR = "results/01_prepare/02_fastq"
OUTDIR = "results/02_demux"
RUN_CELLS = RUN_CELLS

rule all:
    input:
        expand(OUTDIR + "/01_barcode/{run}_bc5.fasta", run=RUNS),
        # expand(OUTDIR + "/02_nanoplexer_5/{run}", run=RUNS),
        # expand(OUTDIR + "/03_nanoplexer_3/{run}", run=RUNS),
        # expand(OUTDIR + "/04_pychopper/{run_cell}.fastq", run_cell=RUN_CELLS),
        expand(OUTDIR + "/05_trim_polya/{run_cell}.fastq.gz", run_cell=RUN_CELLS),

def get_barcode_list_5(run):
    tmp = DAT[DAT["Run"] == run]
    return ["Bar%d" % bc for bc in sorted(set(tmp["Barcode.5"]))]

def get_barcode_list_3(run):
    tmp = DAT[DAT["Run"] == run]
    return ["Bar%d" % bc for bc in sorted(set(tmp["Barcode.3"]))]

rule get_barcode:
    input:
        fa = config["BARCODE_FASTA"]
    output:
        fa5 = OUTDIR + "/01_barcode/{run}_bc5.fasta",
        fa3 = OUTDIR + "/01_barcode/{run}_bc3.fasta"
    params:
        barcode_list_5 = lambda wildcards: get_barcode_list_5(wildcards.run),
        barcode_list_3 = lambda wildcards: get_barcode_list_3(wildcards.run)
    shell:
        """
        samtools faidx {input.fa} {params.barcode_list_5} > {output.fa5}
        samtools faidx {input.fa} {params.barcode_list_3} > {output.fa3}
        """

rule nanoplexer_5:
    input:
        fq = FQDIR + "/{run}.fastq.gz",
        fa = rules.get_barcode.output.fa5
    output:
        out = temp(directory(OUTDIR + "/02_nanoplexer_5/{run}"))
    log:
        OUTDIR + "/02_nanoplexer_5/{run}.log"
    threads:
        24
    shell:
        """
        nanoplexer -t {threads} -b {input.fa} -p {output} {input.fq} &> {log}
        rm {output}/unclassified.fastq
        """

rule nanoplexer_3:
    input:
        fqs = rules.nanoplexer_5.output.out,
        fa = rules.get_barcode.output.fa3
    output:
        out = temp(directory(OUTDIR + "/03_nanoplexer_3/{run}"))
    log:
        OUTDIR + "/03_nanoplexer_3/{run}.log"
    threads:
        24
    shell:
        """
        (
        mkdir {output}
        for fq in {input.fqs}/Bar*.fastq; do
            out={output}/`basename $fq .fastq`
            nanoplexer -t {threads} -b {input.fa} -p $out $fq &> ${{out}}.log
            rm ${{out}}/unclassified.fastq
        done ) &> {log}
        """

# pychopper

def get_cell_barcode_5(wildcards):
    return "Bar%s" % DAT[DAT["Cell"] == wildcards.cell]["Barcode.5"].values[0]

def get_cell_barcode_3(wildcards):
    return "Bar%s" % DAT[DAT["Cell"] == wildcards.cell]["Barcode.3"].values[0]

def get_cell_barcode_sequence_5(wildcards):
    return BARCODES[get_cell_barcode_5(wildcards)][0]

def get_cell_barcode_sequence_3(wildcards):
    return BARCODES[get_cell_barcode_3(wildcards)][0]

rule pychopper:
    input:
        fqs = rules.nanoplexer_3.output,
        txt = config["PRIMER_CONFIG"]
    output:
        txt = OUTDIR + "/04_pychopper/{run}/{cell}.primers.txt",
        fq = temp(OUTDIR + "/04_pychopper/{run}/{cell}.fastq"),
        pdf = OUTDIR + "/04_pychopper/{run}/{cell}.pdf",
        txt2 = OUTDIR + "/04_pychopper/{run}/{cell}.stats.txt"
    log:
        OUTDIR + "/04_pychopper/{run}/{cell}.log"
    params:
        bc_5 = lambda wildcards: get_cell_barcode_5(wildcards),
        bc_3 = lambda wildcards: get_cell_barcode_3(wildcards),
        bc_seq_5 = lambda wildcards: get_cell_barcode_sequence_5(wildcards),
        bc_seq_3 = lambda wildcards: get_cell_barcode_sequence_3(wildcards)
    threads:
        8
    shell:
        """
        echo ">MySSP" > {output.txt}
        echo "{params.bc_seq_5}AAGCAGTGGTATCAACGCAGAGTACATGGG" >> {output.txt}
        echo ">MyVNP" >> {output.txt}
        echo "TCAGACGTGTGCTCTTCCGATC{params.bc_seq_3}" >> {output.txt}
        pychopper -m edlib -r {output.pdf} -b {output.txt} -c {input.txt} -t {threads} \
            -S {output.txt2} {input.fqs}/{params.bc_5}/{params.bc_3}.fastq {output.fq} &> {log}
        """

# remove polyA by custom script

rule trim_polya_and_extract_umi:
    input:
        fq = rules.pychopper.output.fq
    output:
        fq = temp(OUTDIR + "/05_trim_polya/{run}/{cell}.fastq"),
        fq2 = OUTDIR + "/05_trim_polya/{run}/{cell}.fastq.gz"
    log:
        OUTDIR + "/05_trim_polya/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        ./scripts/trim_polya_and_extract_umi.py {input.fq} {output.fq} &> {log}
        pigz -p {threads} -c {output.fq} > {output.fq2}
        """
