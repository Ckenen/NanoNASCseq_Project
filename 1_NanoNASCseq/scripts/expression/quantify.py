#!/usr/bin/env python
import sys
import os
import gzip
import pandas as pd


def main():
    f_sqanti3, f_nascent, f_parental, min_reads, min_tc, outdir = sys.argv[1:]
    min_reads, min_tc = int(min_reads), int(min_tc)
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    # load
    
    if f_sqanti3.endswith(".gz"):
        dat1 = pd.read_csv(gzip.open(f_sqanti3), sep="\t")
    else:
        dat1 = pd.read_csv(f_sqanti3, sep="\t")
    dat1.index = [x[4:] for x in dat1["isoform"]]
    dat1.index.name = "RNA"
    dat1 = dat1[["associated_gene", "associated_transcript", "structural_category", "subcategory"]]

    dat2 = pd.read_csv(f_nascent, sep="\t", header=None)
    dat2.columns = ["Chrom", "Start", "End", "Name", "Reads", "Strand", "Type", "Ts", "TCs", "Detail"]
    dat2.index = dat2["Name"]
    dat2.index.name = "RNA"
    dat2 = dat2[["Reads", "Ts", "TCs"]]

    dat3 = pd.read_csv(f_parental, sep="\t", header=None, index_col=0)
    dat3.columns = ["Parental"]
    dat3.index.name = "RNA"

    dat = pd.concat([dat1, dat2, dat3], sort=False, axis=1)
    dat["Type"] = ["Nascent" if tc >= min_tc else "Exists" for tc in dat["TCs"]]
    dat.index.name = "RNA"
    
    # quantify gene
    
    valid_sc_list = ["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"]
    data_g = dict()
    for gid, sc, reads, tc, parental in dat[["associated_gene", "structural_category", "Reads", "TCs", "Parental"]].values:
        if reads < min_reads:
            continue
        if sc not in valid_sc_list:
            continue
        if gid not in data_g:
            data_g[gid] = [0, 0, 0, 0, 0, 0]
        data_g[gid][0] += 1 # total
        if parental == "P":
            data_g[gid][1] += 1
        elif parental == "M":
            data_g[gid][2] += 1
        if tc >= min_tc: # nascent
            data_g[gid][3] += 1
            if parental == "P":
                data_g[gid][4] += 1
            elif parental == "M":
                data_g[gid][5] += 1
                
    # quantify transcript
    
    data_t = dict()
    for tid, sc, reads, tc, parental in dat[["associated_transcript", "structural_category", "Reads", "TCs", "Parental"]].values:
        if reads < min_reads:
            continue
        if sc != "full-splice_match":
            continue
        if tid not in data_t:
            data_t[tid] = [0, 0, 0, 0, 0, 0]
        data_t[tid][0] += 1 # total
        if parental == "P":
            data_t[tid][1] += 1
        elif parental == "M":
            data_t[tid][2] += 1
        if tc >= min_tc:
            data_t[tid][3] += 1 # nascent
            if parental == "P":
                data_t[tid][4] += 1
            elif parental == "M":
                data_t[tid][5] += 1
    
    # table
    
    columns = ["Total", "Total.P", "Total.M", "Nascent", "Nascent.P", "Nascent.M"]
    gids = []
    rows = []
    for gid, values in data_g.items():
        gids.append(gid)
        rows.append(values)
    df_g = pd.DataFrame(rows, index=pd.Index(gids, name="GeneID"), columns=columns)
            
    tids = []
    rows = []
    for tid, values in data_t.items():
        tids.append(tid)
        rows.append(values)
    df_t = pd.DataFrame(rows, index=pd.Index(tids, name="TranscriptID"), columns=columns)

    # output
    
    dat.to_csv(gzip.open(os.path.join(outdir, "metadata.tsv.gz"), "wt"), sep="\t")
    df_g.to_csv(gzip.open(os.path.join(outdir, "quant_gene.tsv.gz"), "wt"), sep="\t")
    df_t.to_csv(gzip.open(os.path.join(outdir, "quant_transcript.tsv.gz"), "wt"), sep="\t")


if __name__ == "__main__":
    main()