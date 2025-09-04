
# Analysis of NASC-seq

The NASC-seq datasets are from Rickard Sandberg lab (GSE128273) and Tang lab

# Library structures of NASC-seq (Tang lab)

## 20220113 and 20220321

    Forward: 5'-AAGCAGTGGTATCAACGCAGAGTAC-GGG-[RNA]NBAAAAAAAAAAAAAAAAAAAA-GTACTCTGCGTTGATACCACTGCTT-3'
    Reverse: 5'-AAGCAGTGGTATCAACGCAGAGTAC-TTTTTTTTTTTTTTTTTTTTVN[RNA]-CCC-GTACTCTGCGTTGATACCACTGCTT-3'

    5'-TCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-NNNNNNNNNN-CTGTCTCTTATACACATCTGACGCTGCCGACGA-3'
    5'-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG-NNNNNNNNNN-CTGTCTCTTATACACATCTCCGAGCCCACGAGA-3'

## 20220418

    Forward: 5’-AAGCAGTGGTATCAACGCAGAGTAC-ATGGG-[RNA]NBAAAAAAAAAAAAAAAAAAAA-NNNNNNNN-GTACTCTGCGTTGATACCACTGCTT-3’
    Reverse: 5’-AAGCAGTGGTATCAACGCAGAGTAC-NNNNNNNN-TTTTTTTTTTTTTTTTTTTTVN[RNA]-CCCAT-GTACTCTGCGTTGATACCACTGCTT-3’

# GSE128273 NASC-seq (single-end)

    K562 (50uM, 3h), single-end, unstranded


./scripts/stat_category.py -b results/03_mismatch/02_marked_nascent/20220113_NASCseq_K562_04_02_50uM_180min.bam -g ~/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz -t 24 -o analysis/20220113_NASCseq_K562_04_02_50uM_180min.tsv

./scripts/stat_category.py results/03_mismatch/02_marked_nascent/20220113_NASCseq_K562_05_03_50uM_180min.bam ~/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz 24 analysis/20220113_NASCseq_K562_05_03_50uM_180min.tsv

./scripts/stat_category.py ../A01_PROseq_K562/results/06_bam_filtered/GSM2361443_K562_PROseq_NHS_Rep2.bam ~/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz 24 analysis/GSM2361443_K562_PROseq_NHS_Rep2.tsv

./scripts/stat_category.py ../A01_PROseq_K562/results/06_bam_filtered/GSM2545324_K562_PROseq_preseq_0min_plus_rep1.bam ~/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz 24 analysis/GSM2545324_K562_PROseq_preseq_0min_plus_rep1.tsv

./scripts/stat_category.py ../A03_NASCseq_K562_PE/results/03_mismatch/02_marked_nascent/SRR8723979_NASCseq_K562_50uM_60min.bam ~/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz 24 analysis/SRR8723979_NASCseq_K562_50uM_60min.tsv

./scripts/stat_category.py ../A03_NASCseq_K562_PE/results/03_mismatch/02_marked_nascent/SRR8723975_NASCseq_K562_50uM_15min.bam ~/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz 24 analysis/SRR8723975_NASCseq_K562_50uM_15min.tsv


infile1 = "../1_FLAIRseq/results/expression/collapsed/20220719_K562R1/20220719_K562R1.C26.bed.gz"
infile2 = "../1_FLAIRseq/results/mismatch/ratio_consensus/20220719_K562R1/20220719_K562R1.C26.events.tsv"
bedfile = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz"

./scripts/stat_category.py -l ../1_FLAIRseq/results/expression/collapsed/20220719_K562R1/20220719_K562R1.C26.bed.gz -e ../1_FLAIRseq/results/mismatch/ratio_consensus/20220719_K562R1/20220719_K562R1.C26.events.tsv -g /home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz -o a.tsv
