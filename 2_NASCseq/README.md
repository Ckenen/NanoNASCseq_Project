
# Analysis of NASC-seq

The NASC-seq datasets are from Rickard Sandberg lab (GSE128273) and Tanglab

## Workflow

Run the following commands step by step to reproduce the analysis process.

    snakemake -s 1_SnakePrepare.smk -j
    snakemake -s 2_SnakeMapping.smk -j
    snakemake -s 3_SnakeMismatch.smk -j
    snakemake -s 4_SnakeExpression.smk -j
    snakemake -s 5_SnakeStat.smk -j
