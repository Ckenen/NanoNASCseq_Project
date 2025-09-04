#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
FQDIR = "data/datasets"
OUTDIR = "results/2_demux"

rule all:
    input:
        #expand(OUTDIR + "/1_barcodes/{run}.fa", run=RUNS),
        #expand(OUTDIR + "/2_fbilr/{run}.tsv.gz", run=RUNS),
        expand(OUTDIR + "/2_fbilr/{run}.stats.csv.gz", run=RUNS),
        #expand(OUTDIR + "/3_splitted/{run}", run=RUNS),
        #expand(OUTDIR + "/4_trimmed/{run_cell}", run_cell=RUN_CELLS),

rule get_barcodes:
    input:
        fa = config["BARCODES"]
    output:
        fa = OUTDIR + "/1_barcodes/{run}.fa",
        tsv = OUTDIR + "/1_barcodes/{run}.tsv"
    params:
        table = "data/NanoNASCseq.csv"
    run:
        import subprocess
        import pandas as pd
        d = pd.read_csv(params.table)
        d = d[d["Run"] == wildcards.run]
        assert len(d) > 0
        bcs = ["Bar%d" % bc for bc in sorted(set(d["Barcode"]))]
        cmd = "samtools faidx %s %s > %s" % (input.fa, " ".join(bcs), output.fa)
        subprocess.check_call(cmd, shell=True)
        with open(output.tsv, "w+") as fw:
            for cell, bc in d[["Cell", "Barcode"]].values:
                fw.write("%s\tBar%s\n" % (cell, bc))

rule fbilr:
    input:
        fq = FQDIR + "/{run}.fastq.gz",
        fa = rules.get_barcodes.output.fa
    output:
        txt = OUTDIR + "/2_fbilr/{run}.tsv.gz"
    log:
        OUTDIR + "/2_fbilr/{run}.log"
    threads:
        48
    shell:
        """
        ( fbilr -t {threads} -r 10000 -m PE -b {input.fa} {input.fq} | gzip -c > {output.txt} ) &> {log}
        """

rule stat_matrix2:
    input:
        txt = rules.fbilr.output.txt
    output:
        csv = OUTDIR + "/2_fbilr/{run}.stats.csv.gz"
    shell:
        """
        ./scripts/demux/stat_matrix2.py {input.txt} {output.csv}
        """


rule split_reads:
    input:
        fq = FQDIR + "/{run}.fastq.gz",
        mtx = rules.fbilr.output.txt,
        txt = rules.get_barcodes.output.tsv
    output:
        out = directory(OUTDIR + "/3_splitted/{run}")
    log:
        OUTDIR + "/3_splitted/{run}.log"
    threads:
        12
    shell:
        """
        ./scripts/demux/split_reads.py {input} {output} &> {log}
        pigz -p {threads} {output}/*/*.fastq
        """
    
rule trim_reads:
    input:
        fqs = rules.split_reads.output.out
    output:
        out = directory(OUTDIR + "/4_trimmed/{run}/{cell}")
    log:
        OUTDIR + "/4_trimmed/{run}/{cell}.log"
    params:
        fq = rules.split_reads.output.out + "/succeed/{cell}.fastq.gz"
    shell:
        """
        ./scripts/demux/trim_reads.py {params.fq} {output.out} &> {log}
        gzip {output.out}/trimmed.fastq
        """
