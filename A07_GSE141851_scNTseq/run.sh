#!/bin/sh

# output: results/mismatch/04_TC/SRR11683992_K562-4sU-TFEA.tsv_q27.tsv_corrected.tsv
perl ./software/scNT-seq-master/TC_calling/scripts/background_correction.pl \
    -bg results/mismatch/04_TC/SRR11683990_K562-4sU.tsv_q27.tsv \
    -in results/mismatch/04_TC/SRR11683992_K562-4sU-TFEA.tsv_q27.tsv

# output: SRR11683992_K562-4sU-TFEA.TagTC.corrected.bam
perl ./software/scNT-seq-master/TC_calling/scripts/TagIntronicRead_V5.pl \
    -read results/mismatch/04_TC/SRR11683992_K562-4sU-TFEA.tsv_q27.tsv_corrected.tsv \
    -bam results/dropseq/12_DetectedErrors/SRR11683992_K562-4sU-TFEA.bam

# output: SRR11683992_K562-4sU-TFEA.TagTC.corrected_gene_cell_UMI_read.txt
perl ./software/scNT-seq-master/TC_calling/scripts/extract_digital_expression_matrix.pl \
    SRR11683992_K562-4sU-TFEA.TagTC.corrected.bam

# output: SRR11683992_K562-4sU-TFEA_TC_matrix.rds
Rscript ./software/scNT-seq-master/TC_calling/scripts/Generate_T_C_matrix.R \
    SRR11683992_K562-4sU-TFEA.TagTC.corrected_gene_cell_UMI_read.txt \
    600 SRR11683992_K562-4sU-TFEA_TC_matrix.rds




# output: results/mismatch/04_TC/SRR11683995_K562-4sU-TFEA-2ndSS.tsv_q27.tsv_corrected.tsv
perl ./software/scNT-seq-master/TC_calling/scripts/background_correction.pl \
    -bg results/mismatch/04_TC/SRR11683993_K562-4sU-2ndSS.tsv_q27.tsv \
    -in results/mismatch/04_TC/SRR11683995_K562-4sU-TFEA-2ndSS.tsv_q27.tsv

# output: SRR11683995_K562-4sU-TFEA-2ndSS.TagTC.corrected.bam
perl ./software/scNT-seq-master/TC_calling/scripts/TagIntronicRead_V5.pl \
    -read results/mismatch/04_TC/SRR11683995_K562-4sU-TFEA-2ndSS.tsv_q27.tsv_corrected.tsv \
    -bam results/dropseq/12_DetectedErrors/SRR11683995_K562-4sU-TFEA-2ndSS.bam

# output: SRR11683995_K562-4sU-TFEA-2ndSS.TagTC.corrected_gene_cell_UMI_read.txt
perl ./software/scNT-seq-master/TC_calling/scripts/extract_digital_expression_matrix.pl \
    SRR11683995_K562-4sU-TFEA-2ndSS.TagTC.corrected.bam

# output: SRR11683995_K562-4sU-TFEA-2ndSS_TC_matrix.rds
Rscript ./software/scNT-seq-master/TC_calling/scripts/Generate_T_C_matrix.R \
    SRR11683995_K562-4sU-TFEA-2ndSS.TagTC.corrected_gene_cell_UMI_read.txt \
    600 SRR11683995_K562-4sU-TFEA-2ndSS_TC_matrix.rds