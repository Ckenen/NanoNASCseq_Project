#!/usr/bin/env runsnakemake
configfile: "config.yaml"
import glob
SAMPLES = config["SAMPLES"]
FQDIR = "results/01_prepare/03_split"
SAMPLE_BATCH_LIST = []
SAMPLE_TO_BATCH = dict()
for sample in SAMPLES:
    SAMPLE_TO_BATCH[sample] = list()
    d = FQDIR + "/%s" % sample
    if os.path.exists(d):
        for path in sorted(glob.glob(d + "/*.fastq.gz")):
            batch = path.split("/")[-1][:-9]
            SAMPLE_BATCH_LIST.append("%s/%s" % (sample, batch))
            SAMPLE_TO_BATCH[sample].append(batch)
TALLYNN_ROOT = "./software/TallyNN-master"
OUTDIR = "results/02_tallynn"

rule all:
    input:
        expand(OUTDIR + "/01_complement_polyA/{sample_batch}.fastq", sample_batch=SAMPLE_BATCH_LIST),
        expand(OUTDIR + "/02_identify_perfect/{sample_batch}_unambiguous_barcode_R1.fastq", sample_batch=SAMPLE_BATCH_LIST),
        expand(OUTDIR + "/03_identify_perfect_merged/{sample}_unambiguous_barcode_R1.fastq", sample=SAMPLES),
        expand(OUTDIR + "/04_umitools/{sample}.whitelist.txt", sample=SAMPLES),
        expand(OUTDIR + "/05_correct_barcode/{sample_batch}_unambiguous_fixed_barcode_R1.fastq", sample_batch=SAMPLE_BATCH_LIST),
        expand(OUTDIR + "/06_correct_barcode_merged/{sample}_unambiguous_fixed_barcode_R1.fastq", sample=SAMPLES),
        expand(OUTDIR + "/07_merge_full/{sample}_R1.fastq.gz", sample=SAMPLES),
        expand(OUTDIR + "/08_extract_umibc_readname/{sample}.fastq.1.gz", sample=SAMPLES),
        expand(OUTDIR + "/09_umitools_2/{sample}.whitelist.txt", sample=SAMPLES),
        expand(OUTDIR + "/10_umitools_extract/{sample}_R1.fastq.gz", sample=SAMPLES),
        expand(OUTDIR + "/11_remove_adapter/{sample}.fastq.gz", sample=SAMPLES),
        expand(OUTDIR + "/12_reverse_complement/{sample}.fastq.gz", sample=SAMPLES),
        expand(OUTDIR + "/13_remove_polya/{sample}.fastq.gz", sample=SAMPLES),

# Find the polyA in read and reverse complement the read sequences.

rule complement_polyA:
    input:
        fq = FQDIR + "/{sample}/{batch}.fastq.gz"
    output:
        fq = OUTDIR + "/01_complement_polyA/{sample}/{batch}.fastq"
    log:
        OUTDIR + "/01_complement_polyA/{sample}/{batch}.log"
    shell:
        """
        ./scripts/complement_polyA.py --infile {input} --outname {output} &> {log}
        """

# Identify perfect barcode [AABBCCDDEEFFGGHHIIJJKKLL] and imperfect barcode.
# R1 contains barcode and UMI
# R2 contains RNA sequences

rule identify_perfect:
    input:
        fq = rules.complement_polyA.output.fq
    output:
        fq1 = OUTDIR + "/02_identify_perfect/{sample}/{batch}_unambiguous_barcode_R1.fastq",
        fq2 = OUTDIR + "/02_identify_perfect/{sample}/{batch}_unambiguous_barcode_R2.fastq",
        fq3 = OUTDIR + "/02_identify_perfect/{sample}/{batch}_ambiguous_barcode_R1.fastq",
        fq4 = OUTDIR + "/02_identify_perfect/{sample}/{batch}_ambiguous_barcode_R2.fastq",
        txt = OUTDIR + "/02_identify_perfect/{sample}/{batch}.whitelist.txt"
    log:
        OUTDIR + "/02_identify_perfect/{sample}/{batch}.log"
    params:
        prefix = OUTDIR + "/02_identify_perfect/{sample}/{batch}"
    conda:
        "tallynn"
    shell:
        """
        python {TALLYNN_ROOT}/tallynn/python/identify_perfect_nano.py \
            --outname={params.prefix} --infile={input.fq} --whitelist={output.txt} &> {log}
        """

# Merge unambiguous reads

rule merge_unambigous:
    input:
        fqs1 = lambda wildcards: [OUTDIR + "/02_identify_perfect/{sample}/%s_unambiguous_barcode_R1.fastq" % x for x in SAMPLE_TO_BATCH[wildcards.sample]],
        fqs2 = lambda wildcards: [OUTDIR + "/02_identify_perfect/{sample}/%s_unambiguous_barcode_R2.fastq" % x for x in SAMPLE_TO_BATCH[wildcards.sample]]
    output:
        fq1 = OUTDIR + "/03_identify_perfect_merged/{sample}_unambiguous_barcode_R1.fastq",
        fq2 = OUTDIR + "/03_identify_perfect_merged/{sample}_unambiguous_barcode_R2.fastq"
    shell:
        """
        cat {input.fqs1} > {output.fq1}
        cat {input.fqs2} > {output.fq2}
        """

# Generate confident barcode.

rule umitools_whitelist:
    input:
        fq = rules.merge_unambigous.output.fq1
    output:
        txt = OUTDIR + "/04_umitools/{sample}.whitelist.txt"
    log:
        OUTDIR + "/04_umitools/{sample}.whitelist.log"
    params:
        cell_number = lambda wildcards: int(wildcards.sample.split("_")[-1][:-5])
    conda:
        "tallynn"
    shell:
        """
        umi_tools whitelist --stdin={input.fq} --bc-pattern=CCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNN \
            --set-cell-number={params.cell_number} -L {log} > {output.txt}
        """

# Correct ambiguous barcodes to unambiguous barcodes

rule correct_barcode:
    input:
        fq1 = rules.identify_perfect.output.fq3,
        fq2 = rules.identify_perfect.output.fq4,
        txt = rules.umitools_whitelist.output.txt
    output:
        fq1 = OUTDIR + "/05_correct_barcode/{sample}/{batch}_unambiguous_fixed_barcode_R1.fastq",
        fq2 = OUTDIR + "/05_correct_barcode/{sample}/{batch}_unambiguous_fixed_barcode_R2.fastq"
    params:
        prefix = OUTDIR + "/05_correct_barcode/{sample}/{batch}"
    conda:
        "tallynn"
    shell:
        """
        python ./scripts/correct_barcode_nano.py \
            --whitelist={input.txt} --read1={input.fq1} --read2={input.fq2} \
            --outname={params.prefix} --distance=4
        """

# Merge barcode corrected reads

rule merge_corrected:
    input:
        fqs1 = lambda wildcards: [OUTDIR + "/05_correct_barcode/{sample}/%s_unambiguous_fixed_barcode_R1.fastq" % x for x in SAMPLE_TO_BATCH[wildcards.sample]],
        fqs2 = lambda wildcards: [OUTDIR + "/05_correct_barcode/{sample}/%s_unambiguous_fixed_barcode_R2.fastq" % x for x in SAMPLE_TO_BATCH[wildcards.sample]]
    output:
        fq1 = OUTDIR + "/06_correct_barcode_merged/{sample}_unambiguous_fixed_barcode_R1.fastq",
        fq2 = OUTDIR + "/06_correct_barcode_merged/{sample}_unambiguous_fixed_barcode_R2.fastq"
    shell:
        """
        cat {input.fqs1} > {output.fq1}
        cat {input.fqs2} > {output.fq2}
        """

# Merge unambiguous and barcode corrected ambiguous reads.

rule merge_full:
    input:
        fq1 = rules.merge_unambigous.output.fq1,
        fq2 = rules.merge_unambigous.output.fq2,
        fq3 = rules.merge_corrected.output.fq1,
        fq4 = rules.merge_corrected.output.fq2
    output:
        fq1 = OUTDIR + "/07_merge_full/{sample}_R1.fastq.gz",
        fq2 = OUTDIR + "/07_merge_full/{sample}_R2.fastq.gz"
    threads:
        8
    shell:
        """
        cat {input.fq1} {input.fq3} | pigz -p {threads} -c > {output.fq1}
        cat {input.fq2} {input.fq4} | pigz -p {threads} -c > {output.fq2}
        """

# Append the barcode and umi to the read name. (optional)

rule extract_umibc_readname:
    input:
        fq1 = rules.merge_full.output.fq1,
        fq2 = rules.merge_full.output.fq2
    output:
        fq1 = OUTDIR + "/08_extract_umibc_readname/{sample}.fastq.1.gz",
        fq2 = OUTDIR + "/08_extract_umibc_readname/{sample}.fastq.2.gz"
    params:
        prefix = OUTDIR + "/08_extract_umibc_readname/{sample}"
    conda:
        "tallynn"
    shell:
        """
        python {TALLYNN_ROOT}/tallynn/python/extract_umibc_readname.py \
            --read1 {input.fq1} --read2 {input.fq2} --outname {params.prefix}
        """

# Generate barcode whitelist

rule umitools_whitelist_2:
    input:
        fq = rules.merge_full.output.fq1
    output:
        txt = OUTDIR + "/09_umitools_2/{sample}.whitelist.txt"
    log:
        OUTDIR + "/09_umitools_2/{sample}.whitelist.log"
    params:
        cell_number = lambda wildcards: int(wildcards.sample.split("_")[-1][:-5])
    conda:
        "tallynn"
    shell:
        """
        umi_tools whitelist --stdin={input.fq} --bc-pattern=CCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNN \
            --set-cell-number={params.cell_number} -L {log} > {output.txt}
        """

# Extract the barcode and umi from R1 and append to the read name.

rule umitools_extract_bc:
    input:
        fq1 = rules.merge_full.output.fq1,
        fq2 = rules.merge_full.output.fq2,
        txt = rules.umitools_whitelist_2.output.txt
    output:
        fq1 = OUTDIR + "/10_umitools_extract/{sample}_R1.fastq.gz",
        fq2 = OUTDIR + "/10_umitools_extract/{sample}_R2.fastq.gz"
    log:
        OUTDIR + "/10_umitools_extract/{sample}.log"
    conda:
        "tallynn"
    shell:
        """
        umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNN --stdin {input.fq1} \
            --stdout={output.fq1} --read2-in {input.fq2} --read2-out={output.fq2} --whitelist={input.txt} &> {log}
        """

# Remove adapter sequences

rule remove_adapter:
    input:
        fq = rules.umitools_extract_bc.output.fq2
    output:
        fq = OUTDIR + "/11_remove_adapter/{sample}.fastq.gz"
    log:
        OUTDIR + "/11_remove_adapter/{sample}.log"
    threads:
        12
    conda:
        "cutadapt"
    shell:
        """
        cutadapt -j {threads} -e 0.2 -m 200 -a ACTCTGCGTTGATACCACTGCTT -o {output.fq} {input.fq} &> {log}
        """

# Reverse cDNA to RNA sequences.

rule reverse_complement:
    input:
        fq = rules.remove_adapter.output.fq
    output:
        fq = OUTDIR + "/12_reverse_complement/{sample}.fastq.gz"
    threads:
        8
    shell:
        """
        zcat {input.fq} | ./scripts/reverse_fastq.py | pigz -p {threads} -c > {output.fq}
        """

# Remove polyA

rule remove_polya:
    input:
        fq = rules.reverse_complement.output.fq
    output:
        fq = OUTDIR + "/13_remove_polya/{sample}.fastq.gz"
    log:
        OUTDIR + "/13_remove_polya/{sample}.log"
    threads:
        12
    conda:
        "cutadapt"
    shell:
        """
        cutadapt -j {threads} -e 0.2 -m 200 -a AAAAAAAAAAAAAAA -o {output.fq} {input.fq} &> {log}
        """
