#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
GTFDIR = "results/6_assembly/1_stringtie"
OUTDIR = "results/7_assembly_custom"
GROUPS = ['K562', 'mESC', 'MouseBlastocyst']

rule all:
    input:
        expand(OUTDIR + "/1_config/{group}.config.tsv", group=GROUPS),
        expand(OUTDIR + "/2_merged/{group}", group=GROUPS),
        expand(OUTDIR + "/3_sqanti3/{group}", group=GROUPS),
        expand(OUTDIR + "/4_gtf/{group}.all.gtf", group=GROUPS),
        expand(OUTDIR + "/5_gtf_full/{group}.gtf", group=GROUPS),
        expand(OUTDIR + "/5_gtf_full/{group}.bed.gz", group=GROUPS),
        # expand(OUTDIR + "/6_sqanti3_new/{run_cell}", run_cell=run_cells_cellline),

rule make_config:
    output:
        tsv = OUTDIR + "/1_config/{group}.config.tsv"
    run:
        pass
        # cells = get_group_cells(wildcards.group)
        # with open(output.tsv, "w+") as fw:
        #     fw.write("Cell\tGtf\n")
        #     for cell in cells:
        #         path = GTFDIR + "/%s/%s.gtf" % (cell.split(".")[0], cell)
        #         fw.write("%s\t%s\n" % (cell, path))

rule merge_isoforms:
    input:
        tsv = rules.make_config.output.tsv
    output:
        out = directory(OUTDIR + "/2_merged/{group}")
    shell:
        """
        ./scripts/assembly/merge_isoforms.py {input.tsv} {output.out}
        """

rule sqanti3:
    input:
        gtf1 = rules.merge_isoforms.output.out,
        gtf2 = lambda wildcards: get_annotation_gtf(get_group_cells(wildcards.group)[0]),
        fa = lambda wildcards: get_genome_fasta(get_group_cells(wildcards.group)[0])
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
        ./scripts/assembly/run_sqanti3_orf.sh {params.gtf1} {input.gtf2} {input.fa} {threads} {output} &> {log}
        """

rule merge_gtf:
    input:
        sqdir = rules.sqanti3.output.out,
        ref = lambda wildcards: get_annotation_gtf(get_group_cells(wildcards.group)[0]),
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

# def get_group_novel_gtf(cell):
#     group = DAT[DAT["Cell"] == cell]["group"].values[0]
#     if group == "K562":
#         return OUTDIR + "/4_gtf/K562.all.gtf"
#     elif group == "mESC":
#         return OUTDIR + "/4_gtf/mESC.all.gtf"
#     else:
#         return OUTDIR + "/4_gtf/MouseBlastocyst.all.gtf"
    
# rule sqanti3_new:
#     input:
#         gtf1 = GTFDIR + "/{run}/{cell}.gtf",
#         gtf2 = lambda wildcards: get_group_novel_gtf(wildcards.cell),
#         fasta = lambda wildcards: get_genome_fasta(wildcards.cell)
#     output:
#         out = directory(OUTDIR + "/6_sqanti3_new/{run}/{cell}")
#     log:
#         OUTDIR + "/6_sqanti3_new/{run}/{cell}.log"
#     threads:
#         1
#     conda:
#         "SQANTI3.env"
#     shell:
#         """
#         ./scripts/assembly/run_sqanti3_clean.sh {input} {threads} {output} &> {log}
#         """