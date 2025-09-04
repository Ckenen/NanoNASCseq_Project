# Analysis of RNA-seq of K562 under ActD treatment

The cDNA libraries for RNA-seq were stranded (RF).


## Workflow

Run the following commands step by step to reproduce the analysis process.

    snakemake -s 1_SnakePrepare.smk -j
    snakemake -s 2_SnakeMapping.smk -j
    snakemake -s 3_SnakeExpression.smk -j
    snakemake -s 4_SnakeStat.smk -j
    snakemake -s 5_SnakeSNP.smk -j


