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

def get_species(cell):
    return DAT.loc[cell]["Species"]
# def get_group(cell):
#     return DAT.loc[cell]["Group"]
def get_cell_line(cell):
    return DAT.loc[cell]["CellLine"]
def get_label(cell):
    return DAT.loc[cell]["Label"]

DAT_SELECTED = DAT[DAT["Run"].isin(RUNS)]
# DAT_K562 = DAT_SELECTED[DAT_SELECTED["Group"] == "K562"]
# DAT_MESC = DAT_SELECTED[DAT_SELECTED["Group"] == "mESC"]
# DAT_BLASTOCYST = DAT_SELECTED[DAT_SELECTED["Group"] == "MouseBlastocyst"]
RUN_CELLS_SELECTED = list(DAT_SELECTED["RunCell"])
# RUN_CELLS_K562 = list(DAT_K562["RunCell"])
# RUN_CELLS_MESC = list(DAT_MESC["RunCell"])
# RUN_CELLS_BLASTOCYST = list(DAT_BLASTOCYST["RunCell"])
RUN_CELLS = RUN_CELLS_SELECTED
# RUN_CELLS_CELLLINE = RUN_CELLS_K562 + RUN_CELLS_MESC

# GROUPS = ["K562", "mESC", "MouseBlastocyst"]
# TAMA_ROOT = "/home/chenzonggui/software/tama" # Root directory of TAMA

def get_group_cells(group):
    d = DAT
    d = d[(d["Group"] == group) & (d["Time"] == 3) & (d["ActD"].isna())]
    if group == "K562":
        d = d[(d["s4U"] == 0) | (d["s4U"] == 50)]
    elif group == "mESC":
        d = d[(d["s4U"] == 0) | (d["s4U"] == 400)]
    elif group == "MouseBlastocyst":
        d = d[(d["s4U"] == 0) | (d["s4U"] == 400)]
    else:
        assert False
    return list(sorted(d["Cell"]))

def get_estimate_pe_model(cell):
    return "reports/Estimate.Pe.Model.consensus.%s.tsv" % get_cell_line(cell)

def get_chroms(cell):
    species = get_species(cell)
    if species == "Human":
        return ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]
    elif species == "Mouse":
        return ["chr%d" % i for i in range(1, 20)] + ["chrX", "chrY"]

def get_genome_fasta(cell):
    # if "_EXOSC2_AID_" in cell:
    #     return config["EXOSC2AID_GENOME_FASTA"]
    return config["%s_GENOME_FASTA" % get_species(cell).upper()]

def get_genome_splice_mmi(cell):
    # if "_EXOSC2_AID_" in cell:
    #     return config["EXOSC2AID_SPLICE_MMI"]
    return config["%s_SPLICE_MMI" % get_species(cell).upper()]

def get_transcript_bed(cell):
    # if "_EXOSC2_AID_" in cell:
    #     return config["EXOSC2AID_TRANSCRIPT_BED"]
    return config["%s_TRANSCRIPT_BED" % get_species(cell).upper()]

def get_transcript_bed_gz(cell):
    # if "_EXOSC2_AID_" in cell:
    #     return config["EXOSC2AID_TRANSCRIPT_BED_GZ"]
    return config["%s_TRANSCRIPT_BED_GZ" % get_species(cell).upper()]

def get_annotation_gtf(cell):
    # if "_EXOSC2_AID_" in cell:
    #     return config["EXOSC2AID_ANNOTATION_GTF"]
    return config["%s_ANNOTATION_GTF" % get_species(cell).upper()]

def get_transcript_info_tsv(cell):
    return config["%s_TRANSCRIPT_INFO_TSV" % get_species(cell).upper()]

def get_snp_bed(cell):
    # if "_EXOSC2_AID_" in cell:
    #     return config["EXOSC2AID_SNP_BED_GZ"]
    return config["%s_SNP_BED_GZ" % get_species(cell).upper()]

def get_snp_vcf(cell):
    # if "_EXOSC2_AID_" in cell:
    #     assert False
    return config["%s_SNP_VCF_GZ" % get_species(cell).upper()]
