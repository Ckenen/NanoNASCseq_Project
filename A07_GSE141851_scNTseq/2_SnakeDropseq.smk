#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
SAMPLES = config["SAMPLES"]
DROPSEQ_JAR = config["DROPSEQ_JAR"]
SCRIPT_ROOT = config["SCRIPT_ROOT"]
FQDIR = "results/01_prepare/02_fastq"
OUTDIR = "results/02_dropseq"

rule all:
    input:
        # expand(OUTDIR + "/01_FastqToSam/{sample}.bam", sample=SAMPLES),
        # expand(OUTDIR + "/02_TaggedBarcode/{sample}.bam", sample=SAMPLES),
        # expand(OUTDIR + "/03_TaggedMolecular/{sample}.bam", sample=SAMPLES),
        # expand(OUTDIR + "/04_Filtered/{sample}.bam", sample=SAMPLES),
        # expand(OUTDIR + "/05_TrimmedTSO/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/06_TrimmedPolyA/{sample}.bam", sample=SAMPLES),
        # expand(OUTDIR + "/07_SamToFastq/{sample}.fastq", sample=SAMPLES),
        expand(OUTDIR + "/08_StarMapped/{sample}", sample=SAMPLES),
        # expand(OUTDIR + "/09_SortedBam/{sample}.bam", sample=SAMPLES),
        # expand(OUTDIR + "/10_MergedBam/{sample}.bam", sample=SAMPLES),
        expand(OUTDIR + "/11_TaggedGene/{sample}.GeneExonTagged.bam", sample=SAMPLES),
        expand(OUTDIR + "/11_TaggedGene/{sample}.GeneExonTagged.TagIntronic.bam", sample=SAMPLES),
        expand(OUTDIR + "/12_DetectedErrors/{sample}.bam", sample=SAMPLES),

rule FastqToSam: # Sorting by name
    input:
        fq1 = FQDIR + "/{sample}_1.fastq.gz",
        fq2 = FQDIR + "/{sample}_2.fastq.gz",
    output:
        bam = temp(OUTDIR + "/01_FastqToSam/{sample}.bam")
    log:
        OUTDIR + "/01_FastqToSam/{sample}.log"
    conda:
        "picard"
    shell:
        """
        picard FastqToSam \
            FASTQ={input.fq1} \
            FASTQ2={input.fq2} \
            QUALITY_FORMAT=Standard \
            OUTPUT={output.bam} \
            SAMPLE_NAME={wildcards.sample} \
            SORT_ORDER=queryname &> {log}
        """
        
rule ParseCellBarcode:
    input:
        bam = rules.FastqToSam.output.bam
    output:
        bam = temp(OUTDIR + "/02_TaggedBarcode/{sample}.bam"),
        txt = OUTDIR + "/02_TaggedBarcode/{sample}.txt"
    log:
        OUTDIR + "/02_TaggedBarcode/{sample}.log"
    shell:
        """
        java -jar {DROPSEQ_JAR} TagBamWithReadSequenceExtended \
            SUMMARY={output.txt} \
            BASE_RANGE=1-12 \
            BASE_QUALITY=10 \
            BARCODED_READ=1 \
            DISCARD_READ=False \
            TAG_NAME=XC \
            NUM_BASES_BELOW_QUALITY=1 \
            INPUT={input.bam} \
            OUTPUT={output.bam} &> {log}
        """

rule ParseMolecular: # Remove read 1 (DISCARD_READ=True)
    input:
        bam = rules.ParseCellBarcode.output.bam
    output:
        bam = temp(OUTDIR + "/03_TaggedMolecular/{sample}.bam"),
        txt = OUTDIR + "/03_TaggedMolecular/{sample}.txt"
    log:
        OUTDIR + "/03_TaggedMolecular/{sample}.log"
    shell:
        """
        java -jar {DROPSEQ_JAR} TagBamWithReadSequenceExtended \
            SUMMARY={output.txt} \
            BASE_RANGE=13-20 \
            BASE_QUALITY=10 \
            BARCODED_READ=1 \
            DISCARD_READ=True \
            TAG_NAME=XM \
            NUM_BASES_BELOW_QUALITY=1 \
            INPUT={input.bam} \
            OUTPUT={output.bam} &> {log}
        """

rule FilterBam:
    input:
        bam = rules.ParseMolecular.output.bam
    output:
        bam = temp(OUTDIR + "/04_Filtered/{sample}.bam")
    log:
        OUTDIR + "/04_Filtered/{sample}.log"
    shell:
        """
        java -jar {DROPSEQ_JAR} FilterBAM \
            TAG_REJECT=XQ \
            INPUT={input.bam} \
            OUTPUT={output.bam} &> {log}
        """

rule TrimTSO:
    input:
        bam = rules.FilterBam.output.bam
    output:
        bam = temp(OUTDIR + "/05_TrimmedTSO/{sample}.bam"),
        txt = OUTDIR + "/05_TrimmedTSO/{sample}.txt"
    log:
        OUTDIR + "/05_TrimmedTSO/{sample}.log"
    shell:
        """
        java -jar {DROPSEQ_JAR} TrimStartingSequence \
            OUTPUT_SUMMARY={output.txt} \
            SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
            MISMATCHES=0 \
            NUM_BASES=5 \
            INPUT={input.bam} \
            OUTPUT={output.bam} &> {log}
        """

rule TrimPolyA:
    input:
        bam = rules.TrimTSO.output.bam
    output:
        bam = OUTDIR + "/06_TrimmedPolyA/{sample}.bam",
        txt = OUTDIR + "/06_TrimmedPolyA/{sample}.txt"
    log:
        OUTDIR + "/06_TrimmedPolyA/{sample}.log"
    shell:
        """
        java -jar {DROPSEQ_JAR} PolyATrimmer \
            OUTPUT_SUMMARY={output.txt} \
            MISMATCHES=0 \
            NUM_BASES=6 \
            INPUT={input.bam} \
            OUTPUT={output.bam} &> {log}
        """

rule SamToFastq:
    input:
        bam = rules.TrimPolyA.output.bam
    output:
        fq = temp(OUTDIR + "/07_SamToFastq/{sample}.fastq"),
    log:
        OUTDIR + "/07_SamToFastq/{sample}.log"
    conda:
        "picard"
    shell:
        """
        picard SamToFastq INPUT={input.bam} FASTQ={output.fq} &> {log}
        """

rule StarMapping:
    input:
        fq = rules.SamToFastq.output.fq,
        idx = lambda wildcards: get_star_genome(wildcards.sample)
    output:
        directory(OUTDIR + "/08_StarMapped/{sample}")
    log:
        log = OUTDIR + "/08_StarMapped/{sample}.log"
    threads:
        12
    conda:
        "star"
    shell:
        """
        mkdir -p {output}
        STAR --runThreadN {threads} \
            --genomeDir {input.idx} \
            --outFilterScoreMinOverLread 0.3 \
            --outFilterMatchNminOverLread 0.3 \
            --readFilesCommand zcat \
            --readFilesIn {input.fq} \
            --outFileNamePrefix {output}/{wildcards.sample}. &> {log}
        """

rule SortSam: # sorted by name
    input:
        samdir = rules.StarMapping.output
    output:
        bam = temp(OUTDIR + "/09_SortedBam/{sample}.bam")
    log:
        OUTDIR + "/09_SortedBam/{sample}.log"
    conda:
        "picard"
    shell:
        """
        picard SortSam \
            INPUT={input.samdir}/{wildcards.sample}.Aligned.out.sam \
            OUTPUT={output.bam} \
            SORT_ORDER=queryname &> {log}
        """

rule MergeBam: # Added MD tag, sorted by coordinates
    input:
        bam1 = rules.TrimPolyA.output.bam,
        bam2 = rules.SortSam.output.bam,
        fasta = lambda wildcards: get_genome_fasta(wildcards.sample)
    output:
        bam = OUTDIR + "/10_MergedBam/{sample}.bam"
    log:
        OUTDIR + "/10_MergedBam/{sample}.log"
    conda:
        "picard"
    shell:
        """
        picard MergeBamAlignment \
            REFERENCE_SEQUENCE={input.fasta} \
            UNMAPPED_BAM={input.bam1} \
            ALIGNED_BAM={input.bam2} \
            OUTPUT={output.bam} \
            INCLUDE_SECONDARY_ALIGNMENTS=false \
            PAIRED_RUN=false &> {log}
        """

rule TagExonicRead:
    input:
        bam = rules.MergeBam.output.bam,
        gtf = lambda wildcards: get_annotation_gtf(wildcards.sample)
    output:
        bam = temp(OUTDIR + "/11_TaggedGene/{sample}.GeneExonTagged.bam")
    log:
        OUTDIR + "/11_TaggedGene/{sample}.GeneExonTagged.log"
    conda:
        "picard"
    shell:
        """
        java -jar {DROPSEQ_JAR} TagReadWithGeneExon \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            ANNOTATIONS_FILE={input.gtf} TAG=GE &> {log}
        """

rule TagIntronicRead:
    input:
        bam = rules.TagExonicRead.output.bam,
        gtf = lambda wildcards: get_annotation_gtf(wildcards.sample)
    output:
        bam1 = OUTDIR + "/11_TaggedGene/{sample}.GeneExonTagged.TagIntronic.bam",
        bam2 = OUTDIR + "/11_TaggedGene/{sample}.GeneExonTagged.TagIntronicOnly.bam"
    log:
        OUTDIR + "/11_TaggedGene/{sample}.GeneExonTagged.TagIntronic.log"
    shell:
        """
        perl {SCRIPT_ROOT}/TagIntronicRead_V3.pl \
            -gtf {input.gtf} \
            -bam {input.bam} &> {log} 
        """

rule DetectBeadSynthesisErrors:
    input:
        bam = rules.TagIntronicRead.output.bam1
    output:
        bam = OUTDIR + "/12_DetectedErrors/{sample}.bam",
        txt = OUTDIR + "/12_DetectedErrors/{sample}.stats.txt",
        txt2 = OUTDIR + "/12_DetectedErrors/{sample}.summary.txt",
    log:
        OUTDIR + "/12_DetectedErrors/{sample}.log"
    params:
        num_cells = 1200
    shell:
        """
        java -jar {DROPSEQ_JAR} DetectBeadSynthesisErrors \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            OUTPUT_STATS={output.txt} \
            SUMMARY={output.txt2} \
            NUM_BARCODES={params.num_cells} \
            PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC &> {log}
        """
