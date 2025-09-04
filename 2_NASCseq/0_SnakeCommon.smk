#!/usr/bin/env runsnakemake
configfile: "config.yaml"
RUNS = config["RUNS"]
import pandas as pd
DAT = pd.read_excel("data/NASCseq.xlsx")
DAT = DAT[DAT["Run"].isin(RUNS)]
DAT["RunCell"] = ["%s/%s" % (run, cell) for run, cell in DAT[["Run", "Cell"]].values]
# print(DAT)
RUN_CELLS = DAT["RunCell"]