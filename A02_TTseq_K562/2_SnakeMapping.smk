#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config['SAMPLES']
FQDIR = "results/01_prepare/03_bowtie2_mapping_ribo"
OUTDIR = "results/02_mapping"

rule all:
    input:
        # OUTDIR + "/01_star_index_genome/genome_star_index",
        # expand(OUTDIR + "/02_star_mapping/{sample}", sample=SAMPLES),
        expand(OUTDIR + "/03_bam_filtered/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/03_bam_filtered/{sample}.flagstat", sample=SAMPLES),
        expand(OUTDIR + "/04_markdup/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/04_markdup/{sample}.flagstat", sample=SAMPLES),

rule star_build_genome:
    input:
        fa = config["GENOME_FASTA"],
        gtf = config["GENOME_GTF"]
    output:
        idx = directory(OUTDIR + "/01_star_index_genome/genome_star_index")
    log:
        OUTDIR + "/01_star_index_genome/genome_star_index.log"
    threads:
        12
    conda:
        "star"
    shell:
        """
        mkdir -p {output.idx}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.idx} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} &> {log}
        """

rule star_mapping:
    input:
        fq1 = FQDIR + "/{sample}.unmapped.fastq.1.gz",
        fq2 = FQDIR + "/{sample}.unmapped.fastq.2.gz",
        idx = rules.star_build_genome.output.idx
    output:
        out = directory(OUTDIR + "/02_star_mapping/{sample}")
    conda:
        "star"
    log:
        OUTDIR + "/02_star_mapping/{sample}.log"
    params:
        prefix = OUTDIR + "/02_star_mapping/{sample}/{sample}"
    threads:
        24
    shell:
        """(
        mkdir -p {output.out}
        STAR --runThreadN {threads} \
            --outFileNamePrefix {params.prefix}. \
            --genomeDir {input.idx} \
            --readFilesCommand zcat \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fq1} {input.fq2} ) &> {log}
        """

rule filter_bam:
    input:
        bamdir = rules.star_mapping.output.out
    output:
        bam = OUTDIR + "/03_bam_filtered/{sample}.bam"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} -q 30 -d "NH:1" --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            -f 2 -F 2308 -o {output.bam} {input}/{wildcards.sample}.Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam}
        """

rule markdup:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/04_markdup/{sample}.bam",
        tmpdir = temp(directory(OUTDIR + "/04_markdup/{sample}.TMP"))
    log:
        OUTDIR + "/04_markdup/{sample}.log"
    threads:
        4
    shell:
        """
        mkdir -p {output.tmpdir}
        sambamba markdup -r -t {threads} --tmpdir={output.tmpdir} {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
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