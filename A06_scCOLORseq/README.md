# Analysis of scCOLOR-seq datasets

This repository contains the `Snakemake` workflow to analyze the scCOLOR-seq datasets and statistics of the Nanopore-based single-cell RNA sequencing technology yield.

## Datasets

The 6 datasets were download from GEO under accession GSE162053 (High throughput error correction using dual nucleotide dimer blocks allows direct single-cell nanopore transcriptome sequencing, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162053). 

The information of the datasets were show below:

| Name | Sample | Run | AvgSpotLen | Instrument | Species | Description | 
| :- | :- | :- | :- | :- | :- | :- |
| SRR13120810_3T3_HEK293_1000cells | GSM4932357 | SRR13120810 | 1044 | MinION | Homo sapiens and Mus musculus | 3T3/HEK (1000 cells) |
| SRR13397411_3T3_HEK293_500cells | GSM5012389 | SRR13397411 | 970 | MinION | Homo sapiens and Mus musculus | 3T3/HEK (500 cells) |
| SRR13397412_Myeloma_500cells | GSM5012390 | SRR13397412 | 869 | MinION | Homo sapiens | Myeloma cell lines (500 cells) |
| SRR13397413_Myeloma_500cells | GSM5012391 | SRR13397413 | 1014 | MinION | Homo sapiens | Myeloma cell lines (500 cells) |
| SRR13397414_Ewings_500cells | GSM5012392 | SRR13397414 | 1145 | MinION | Homo sapiens | Ewings cells (500 cells) |
| SRR13509039_Myeloma_1200cells | GSM5031743 | SRR13509039 | 1148 | PromethION | Homo sapiens | Myeloma cell lines (1200 cells) |

## Pipelines

The pipeline includes 3 steps described as follows:

 - Step 1: Download SRA files, converted to FASTQ files, and split the FASTQ for parrallel running.

   snakemake -s 1_SnakePrepare.smk -j --rerun-triggers mtime --use-conda --conda-not-block-search-path-envvars

 - Step 2: Identify barcode and UMI, and remove exogenous sequences. We adopt `pipeline_nanopore.py` in `TallyNN` (https://github.com/Acribbs/TallyNN) into `Snakemake`.

   snakemake -s 2_SnakeTallyNN.smk -j --rerun-triggers mtime --use-conda --conda-not-block-search-path-envvars

 - Step 3: Mapping the clean reads to integrated genome (Homo sapiens and Mus musculus).

   snakemake -s 3_SnakeMapping.smk -j --rerun-triggers mtime --use-conda --conda-not-block-search-path-envvars

## Results

Yield of sequencing data

![Figure](figures/scCOLORseq_recovered_of_reads.png)

Summary:

    Mean:
    TotalRatio              100.000000
    PolyARatio               43.255731
    AdapterRatio             29.875306
    MergedFullRatio          13.132245
    UmitoolsExtractRatio     10.115523
    MappedRatio               7.630489
    FilteredRatio             7.345930

    Std:
    TotalRatio               0.000000
    PolyARatio              10.214190
    AdapterRatio             6.549523
    MergedFullRatio          3.147468
    UmitoolsExtractRatio     3.271037
    MappedRatio              2.746715
    FilteredRatio            2.644300

## References

Philpott, M., Watson, J., Thakurta, A., Brown, T., Jr., Brown, T., Sr., Oppermann, U., & Cribbs, A. P. (2021). Nanopore sequencing of single-cell transcriptomes with scCOLOR-seq. Nat Biotechnol, 39(12), 1517-1520. https://doi.org/10.1038/s41587-021-00965-w 
