#!/usr/bin/env python
import sys
import pandas as pd

d = pd.read_csv(sys.argv[1], sep="\t")
d = d[d["Canonical"]]
d1 = d[d["GeneType"] == "protein_coding"]
d2 = d1[d1["FPKM"] >= 1]
d3 = d1[d1["FPKM.New"] >= 1]
protein_coding = len(set(d2["GeneName"]))
protein_coding_nascent = len(set(d3["GeneName"]))
d1 = d[d["GeneType"] == "lncRNA"]
d2 = d1[d1["FPKM"] >= 1]
d3 = d1[d1["FPKM.New"] >= 1]
lncRNA = len(set(d2["GeneName"]))
lncRNA_nascent = len(set(d3["GeneName"]))
print("Protein_coding\tProtein_coding.Nascent\tlncRNA\tlncRNA.Nascent")
print("\t".join(map(str, [protein_coding, protein_coding_nascent, lncRNA, lncRNA_nascent])))