#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
GTFDIR = "results/6_assembly/1_stringtie"
OUTDIR = "results/7_assembly_novel"
GROUPS = ['K562', 'mESC', 'MouseBlastocyst']

rule all:
    input:
        expand(OUTDIR + "/1_config/{group}.tsv", group=GROUPS),
        expand(OUTDIR + "/2_merged_isoforms/{group}", group=GROUPS),
        expand(OUTDIR + "/3_sqanti3/{group}", group=GROUPS),
        expand(OUTDIR + "/4_gtf/{group}.all.gtf", group=GROUPS),
        expand(OUTDIR + "/5_gtf_full/{group}.bed.gz", group=GROUPS),

def get_group_cells(group):
    df = DAT
    if group == "K562":
        d = df[(df["CellType"] == "K562") & (df["s4U"].isin([0, 50])) & (df["Time"] == 3) & (df["ActD"].isna())]
    elif group == "mESC":
        d = df[(df["CellType"] == "mESC") & (df["s4U"].isin([0, 400])) & (df["Time"] == 3) & (df["ActD"].isna())]
    elif group == "MouseBlastocyst":
        d = df[(df["CellType"] == "MouseBlastocyst") & (df["s4U"].isin([0, 400])) & (df["Time"] == 3) & (df["ActD"].isna())]
    else:
        assert False
    return list(sorted(df["Cell"]))

def get_group_gtfs(group):
    paths = []
    for cell in get_group_cells(group):
        run = cell.split(".")[0]
        paths.append(GTFDIR + "/%s/%s.gtf" % (run, cell))
    return paths

rule make_config:
    output:
        tsv = OUTDIR + "/1_config/{group}.tsv"
    run:
        with open(output.tsv, "w") as fout:
            for cell in get_group_cells(wildcards.group):
                fout.write("%s\t%s\n" % (cell, GTFDIR + "/%s/%s.gtf" % (cell.split(".")[0], cell)))

rule merge_isoforms:
    input:
        gtfs = lambda wildcards: get_group_gtfs(wildcards.group),
        tsv = rules.make_config.output.tsv
    output:
        out = directory(OUTDIR + "/2_merged_isoforms/{group}")
    shell:
        """
        ./scripts/assembly/merge_isoforms.py {input.tsv} {output.out}
        """

rule sqanti3:
    input:
        gtf1 = rules.merge_isoforms.output.out,
        gtf2 = lambda wildcards: config["%s_ANNOTATION_GTF" % ("HUMAN" if wildcards.group == "K562" else "MOUSE")],
        fa = lambda wildcards: config["%s_GENOME_FASTA" % ("HUMAN" if wildcards.group == "K562" else "MOUSE")]
    output:
        out = directory(OUTDIR + "/3_sqanti3/{group}")
    log:
        OUTDIR + "/3_sqanti3/{group}.log"
    params:
        gtf1 = rules.merge_isoforms.output.out + "/merged.gtf"
    threads:
        8
    conda:
        "SQANTI3.env"
    shell:
        """
        ../share/scripts/run_sqanti3_orf.sh {params.gtf1} {input.gtf2} {input.fa} {threads} {output} &> {log}
        """

rule merge_gtf:
    input:
        sqdir = rules.sqanti3.output.out,
        ref = lambda wildcards: config["%s_ANNOTATION_GTF" % ("HUMAN" if wildcards.group == "K562" else "MOUSE")],
    output:
        gtf1 = OUTDIR + "/4_gtf/{group}.novel.gtf",
        gtf2 = OUTDIR + "/4_gtf/{group}.all.gtf"
    params:
        tsv = rules.sqanti3.output.out + "/merged_classification.txt",
        gtf = rules.sqanti3.output.out + "/merged_corrected.gtf.cds.gff",
    shell:
        """
        ./scripts/assembly/report_novel_isoforms.py {params.tsv} {params.gtf} | sort -k1,1 -k4,4n > {output.gtf1}
        cat {input.ref} {output.gtf1} | sort -k1,1 -k4,4n > {output.gtf2}
        bgzip -c {output.gtf1} > {output.gtf1}.gz
        tabix -p gff {output.gtf1}.gz
        bgzip -c {output.gtf2} > {output.gtf2}.gz
        tabix -p gff {output.gtf2}.gz
        """

rule add_attributes:
    input:
        gtf = rules.merge_gtf.output.gtf2
    output:
        gtf = OUTDIR + "/5_gtf_full/{group}.gtf"
    shell:
        """
        ./scripts/assembly/add_attributes.py {input.gtf} {output.gtf}
        bgzip -c {output.gtf} > {output.gtf}.gz
        tabix -p gff {output.gtf}.gz
        """

rule fetch_transcripts:
    input:
        gtf = rules.add_attributes.output.gtf
    output:
        bed = OUTDIR + "/5_gtf_full/{group}.bed.gz",
        tbi = OUTDIR + "/5_gtf_full/{group}.bed.gz.tbi"
    shell:
        """
        fetch_transcripts.py {input.gtf} | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """
