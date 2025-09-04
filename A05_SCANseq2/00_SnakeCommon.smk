#!/usr/bin/env runsnakemake
configfile: "config.yaml"

import pandas as pd
DAT = pd.read_excel("data/SCANseq2.xlsx")
RUN_CELLS = []
for run, cell in DAT[["Run", "Cell"]].values:
    RUN_CELLS.append("%s/%s" % (run, cell))
print("Cells:", len(RUN_CELLS))

BARCODES = dict()
from Bio import SeqIO
for read in SeqIO.parse(config["BARCODE_FASTA"], "fasta"):
    seq1 = str(read.seq)
    seq2 = str(read.seq.reverse_complement())
    BARCODES[read.name] = [seq1, seq2]
# print(barcodes)

def get_species(cell):
    if "_Human_" in cell:
        return "Human"
    elif "_Mouse_" in cell:
        return "Mouse"
    assert False
