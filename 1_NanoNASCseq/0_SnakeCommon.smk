#!/usr/bin/env runsnakemake
configfile: "config.yaml"
import os
import numpy as np
import pandas as pd

RUNS = []
RUNS.extend(config["RUNS_K562"])
RUNS.extend(config["RUNS_K562_ACTD"])
RUNS.extend(config["RUNS_MESC"])
RUNS.extend(config["RUNS_MOUSE_BLASTOCYST"])
RUNS.extend(config["RUNS_MIX"])
RUNS.extend(config["RUNS_K562_FUCCI"])
RUNS.extend(config["RUNS_MESC_EXOSC2"])

THREADS = 8

DAT = pd.read_csv(config["TABLE"])
DAT.index = DAT["Cell"]
DAT["RunCell"] = ["/".join(vs) for vs in DAT[["Run", "Cell"]].values]
RUN_CELLS = DAT['RunCell']

DAT_SELECTED = DAT[DAT["Run"].isin(RUNS)]
RUN_CELLS_SELECTED = list(DAT_SELECTED["RunCell"])
RUN_CELLS = RUN_CELLS_SELECTED

def get_species(cell):
    return DAT.loc[cell]["Species"]

def get_cell_type(cell):
    return DAT.loc[cell]["CellType"]

def get_estimate_pe_model(cell):
    return "reports/Estimate.Pe.Model.consensus.%s.tsv" % get_cell_line(cell)

# def get_chroms(cell):
#     species = get_species(cell)
#     if species == "Human":
#         return ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]
#     elif species == "Mouse":
#         return ["chr%d" % i for i in range(1, 20)] + ["chrX", "chrY"]

