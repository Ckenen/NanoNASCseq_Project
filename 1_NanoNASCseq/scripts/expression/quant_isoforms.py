#!/usr/bin/env python
import sys
from collections import defaultdict
import pandas as pd


def counting_alleles(d):
    counter = defaultdict(int)
    for v in d["Alleles"]:
        if isinstance(v, str):
            for a in v.split(";"):
                counter[a] += 1
    return ";".join(["%s:%d" % (k, v) for k, v in sorted(counter.items())])


def main():
    stat_tsv, event_tsv, allele_tsv, min_read, min_tc, out_tsv = sys.argv[1:]
    min_read = int(min_read)
    min_tc = int(min_tc)

    dat1 = pd.read_csv(stat_tsv, sep="\t")
    selected_index_list = []
    for name, tmp in dat1.groupby("Name"):
        if len(tmp) == 1:
            selected_index_list.append(tmp.index.values[0])
        else:
            cg = list(set(tmp["Category"]))
            if len(cg) == 1 and cg[0] == "FSM":
                tmp = tmp.sort_values(by="TranscriptID")
                selected_index_list.append(tmp.index.values[0])
            
    dat1 = dat1.loc[selected_index_list]
    dat1.index = dat1["Name"]

    dat2 = pd.read_csv(event_tsv, sep="\t", index_col=0)
    dat2 = dat2[["Size", "T-C"]]

    dat3 = pd.read_csv(allele_tsv, sep="\t", index_col=0)
    dat3 = dat3[["Alleles"]]
    
    dat = pd.concat([dat1, dat2, dat3], axis=1, join="inner")
    dat = dat[dat["Size"] >= min_read]
    rows = []
    for tid, tmp in dat.groupby(by="TranscriptID"):
        tmp1 = tmp[tmp["T-C"] >= min_tc]
        row = [tid, len(tmp), len(tmp1), counting_alleles(tmp), counting_alleles(tmp1)]
        rows.append(row)
    m = pd.DataFrame(rows, columns=["TranscriptID", "Total", "New", "Total.Alleles", "New.Alleles"])
    m = m.sort_values(by="TranscriptID")
    m.to_csv(out_tsv, sep="\t", index=False)
    

if __name__ == "__main__":
    main()