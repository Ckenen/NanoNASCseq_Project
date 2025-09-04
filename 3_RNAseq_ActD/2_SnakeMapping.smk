#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
SPECIES = config["SPECIES"]
FQDIR = "results/1_prepare/2_bowtie2"
OUTDIR = "results/2_mapping"

rule all:
    input:
        # OUTDIR + "/1_merged_genome/human_fly.fa",
        # OUTDIR + "/1_merged_genome/human_fly_star_index",
        # expand(OUTDIR + "/2_star/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/3_filtered/{sample}.{species}.bam", sample=SAMPLES, species=SPECIES),
        expand(OUTDIR + "/3_filtered/{sample}.{species}.flagstat", sample=SAMPLES, species=SPECIES),
        expand(OUTDIR + "/4_rmdup/{sample}.{species}.bam", sample=SAMPLES, species=SPECIES),
        expand(OUTDIR + "/4_rmdup/{sample}.{species}.flagstat", sample=SAMPLES, species=SPECIES),

rule merge_genome:
    input:
        fa1 = config["HUMAN_FASTA"],
        fa2 = config["FLY_FASTA"],
        gtf1 = config["HUMAN_GTF"],
        gtf2 = config["FLY_GTF"]
    output:
        fa = OUTDIR + "/1_merged_genome/human_fly.fa",
        fai = OUTDIR + "/1_merged_genome/human_fly.fa.fai",
        gtf = OUTDIR + "/1_merged_genome/human_fly.gtf",
        gtf2 = OUTDIR + "/1_merged_genome/human_fly.gtf.gz",
        tbi2 = OUTDIR + "/1_merged_genome/human_fly.gtf.gz.tbi"
    log:
       OUTDIR + "/1_merged_genome/human_fly.fa.log" 
    shell:
        """
        cat {input.fa1} {input.fa2} > {output.fa}
        samtools faidx {output.fa}
        cat {input.gtf1} {input.gtf2} | grep -v '#' | sort -k1,1 -k4,4n -k5,5n > {output.gtf} 
        bgzip -c {output.gtf} > {output.gtf2}
        tabix -p gff {output.gtf2}
        """

rule star_build_genome:
    input:
        fa = rules.merge_genome.output.fa,
        gtf = rules.merge_genome.output.gtf
    output:
        idx = directory(OUTDIR + "/1_merged_genome/human_fly_star_index")
    log:
        OUTDIR + "/1_merged_genome/human_fly_star_index.log"
    conda:
        "star"
    threads:
        20
    shell:
        """(
        mkdir -p {output}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} ) &> {log}
        """

rule star_mapping:
    input:
        fq1 = FQDIR + "/{sample}.unmapped.fastq.1.gz",
        fq2 = FQDIR + "/{sample}.unmapped.fastq.2.gz",
        idx = rules.star_build_genome.output.idx
    output:
        directory(OUTDIR + "/2_star/{sample}")
    log:
        OUTDIR + "/2_star/{sample}.log"
    conda:
        "star"
    threads:
        20
    shell:
        """(
        mkdir -p {output}
        STAR --runThreadN {threads} \
            --outFileNamePrefix {output}/{wildcards.sample}. \
            --genomeDir {input.idx} \
            --readFilesCommand zcat \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fq1} {input.fq2} ) &> {log}
        """

rule filter_and_split: # filter and split
    input:
        rules.star_mapping.output
    output:
        bam = OUTDIR + "/3_filtered/{sample}.{species}.bam"
    log:
        OUTDIR + "/3_filtered/{sample}.{species}.log"
    params:
        pattern = lambda wildcards: config["HUMAN_SEQNAME_PATTERN"] if wildcards.species == "human" else config["FLY_SEQNAME_PATTERN"]
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} -q 30 -d "NH:1" -f 2 -F 2308 \
            --expr 'rname =~ "{params.pattern}"' -o {output.bam} \
            {input}/{wildcards.sample}.Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule rmdup:
    input:
        bam = rules.filter_and_split.output.bam
    output:
        bam = OUTDIR + "/4_rmdup/{sample}.{species}.bam"
    log:
        OUTDIR + "/4_rmdup/{sample}.{species}.log"
    threads:
        4
    shell:
        """
        sambamba markdup -r -t {threads} {input.bam} {output.bam} &> {log}
        """

# Common rules

rule bam_flagstat:
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