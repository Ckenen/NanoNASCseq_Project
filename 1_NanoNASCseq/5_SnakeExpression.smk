#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
BAMDIR = "results/3_mapping/5_markdup"
EVENTDIR = "results/4_mismatch/2_mismatch_ratios/3_consensus"
OUTDIR = "results/5_expression"
# RUN_CELLS = RUN_CELLS[:1]
CUTOFFS = ["min_read_1_min_tc_1", "min_read_2_min_tc_1", "min_read_2_min_tc_2"]

rule all:
    input:
        expand(OUTDIR + "/1_expressed_alleles/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/2_collapsed/{run_cell}.gtf", run_cell=RUN_CELLS),
        expand(OUTDIR + "/3_sqanti3/{run_cell}", run_cell=RUN_CELLS),
        expand(OUTDIR + "/4_quant_genes/{cutoff}/{run_cell}.tsv", cutoff=CUTOFFS, run_cell=RUN_CELLS),
        expand(OUTDIR + "/5_isoform_category/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/6_quant_isoforms/{cutoff}/{run_cell}.tsv", cutoff=CUTOFFS, run_cell=RUN_CELLS)

rule stat_expressed_alleles:
    input:
        bam = BAMDIR + "/{run}/{cell}.bam",
        vcf = lambda wildcards: config["%s_SNP_VCF_GZ" % get_species(wildcards.cell).upper()]
    output:
        tsv = OUTDIR + "/1_expressed_alleles/{run}/{cell}.tsv"
    log:
        OUTDIR + "/1_expressed_alleles/{run}/{cell}.log"
    shell:
        """
        ./scripts/expression/stat_expressed_alleles.py {input} > {output.tsv} 2> {log}
        """

rule collapse:
    input:
        bam = BAMDIR + "/{run}/{cell}.bam",
    output:
        bed = OUTDIR + "/2_collapsed/{run}/{cell}.bed",
        bed_gz = OUTDIR + "/2_collapsed/{run}/{cell}.bed.gz",
        gtf = OUTDIR + "/2_collapsed/{run}/{cell}.gtf",
        gtf_gz = OUTDIR + "/2_collapsed/{run}/{cell}.gtf.gz"
    log:
        OUTDIR + "/2_collapsed/{run}/{cell}.log"
    shell:
        """(
        ./scripts/expression/collapse_umi.py {input.bam} | sort -k1,1 -k2,2n -k3,3n > {output.bed}
        bgzip -c {output.bed} > {output.bed_gz}
        tabix -p bed -f {output.bed_gz}
        ./scripts/expression/bed2gtf.py {output.bed} | sort -k1,1 -k4,4n > {output.gtf}
        bgzip -c {output.gtf} > {output.gtf_gz}
        tabix -p gff -f {output.gtf_gz} ) &> {log}
        """

rule sqanti3:
    input:
        gtf1 = rules.collapse.output.gtf,
        gtf2 = lambda wildcards: config["%s_ANNOTATION_GTF" % get_species(wildcards.cell).upper()],
        fa = lambda wildcards: config["%s_GENOME_FASTA" % get_species(wildcards.cell).upper()]
    output:
        out = directory(OUTDIR + "/3_sqanti3/{run}/{cell}")
    log:
        OUTDIR + "/3_sqanti3/{run}/{cell}.log"
    conda:
        "SQANTI3.env"
    threads:
        4
    shell:
        """
        ../share/scripts/run_sqanti3_clean.sh {input} {threads} {output} &> {log}
        """

rule quant_genes:
    input:
        sqanti3_out = rules.sqanti3.output.out,
        event_tsv = EVENTDIR + "/{run}/{cell}.events.tsv",
        allele_tsv = rules.stat_expressed_alleles.output.tsv
    output:
        tsv = OUTDIR + "/4_quant_genes/min_read_{size}_min_tc_{tc}/{run}/{cell}.tsv"
    log:
        OUTDIR + "/4_quant_genes/min_read_{size}_min_tc_{tc}/{run}/{cell}.log"
    params:
        sqanti3_tsv = rules.sqanti3.output.out + "/{cell}_classification.txt"
    shell:
        """
        if [ -f "{params.sqanti3_tsv}" ]; then
            ./scripts/expression/quant_genes.py {params.sqanti3_tsv} {input.event_tsv} \
                {input.allele_tsv} {wildcards.size} {wildcards.tc} {output.tsv} &> {log}
        else
            echo -e 'GeneID\\tTotal\\tNew\\tTotal.Alleles\\tNew.Alleles' > {output.tsv}
        fi
        """

rule stat_isoform_category:
    input:
        gtf = lambda wildcards: config["%s_ANNOTATION_GTF" % get_species(wildcards.cell).upper()],
        bed = rules.collapse.output.bed_gz
    output:
        tsv = OUTDIR + "/5_isoform_category/{run}/{cell}.tsv"
    log:
        OUTDIR + "/5_isoform_category/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        ./scripts/expression/stat_isoform_category.py {input} {output} &> {log}
        """

rule quant_isoforms:
    input:
        stat_tsv = rules.stat_isoform_category.output.tsv,
        event_tsv = EVENTDIR + "/{run}/{cell}.events.tsv",
        allele_tsv = rules.stat_expressed_alleles.output.tsv
    output:
        tsv = OUTDIR + "/6_quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.tsv"
    log:
        OUTDIR + "/6_quant_isoforms/min_read_{size}_min_tc_{tc}/{run}/{cell}.log"
    shell:
        """
        ./scripts/expression/quant_isoforms.py {input.stat_tsv} {input.event_tsv} \
            {input.allele_tsv} {wildcards.size} {wildcards.tc} {output.tsv} &> {log}
        """

