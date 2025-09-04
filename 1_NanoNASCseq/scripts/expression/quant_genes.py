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
    sqanti3_tsv, event_tsv, allele_tsv, min_read, min_tc, out_tsv = sys.argv[1:]
    
    min_read = int(min_read)
    min_tc = int(min_tc)
    dat1 = pd.read_csv(sqanti3_tsv, sep="\t", index_col=0)
    dat2 = pd.read_csv(event_tsv, sep="\t", index_col=0)
    dat2 = dat2[["Size", "T-C"]]
    dat3 = pd.read_csv(allele_tsv, sep="\t", index_col=0)
    dat3 = dat3[["Alleles"]]
    # assert len(dat1) == len(dat2)
    # assert len(dat2) == len(dat3)
    
    dat = pd.concat([dat1, dat2, dat3], axis=1, join="inner")
    dat = dat[dat["Size"] >= min_read]
    sc_list = ['full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog']
    dat = dat[[sc in sc_list for sc in dat["structural_category"]]]
    if len(dat) > 0:
        rows = []
        for gid, tmp in dat.groupby(by="associated_gene"):
            tmp1 = tmp[tmp["T-C"] >= min_tc]
            row = [gid, len(tmp), len(tmp1), counting_alleles(tmp), counting_alleles(tmp1)]
            rows.append(row)
        m = pd.DataFrame(rows, columns=["GeneID", "Total", "New", "Total.Alleles", "New.Alleles"])
        m = m.sort_values(by="GeneID")
        m.to_csv(out_tsv, sep="\t", index=False)
    else:
        with open(out_tsv, "w+") as fw:
            fw.write("GeneID\tTotal\tNew\tTotal.Alleles\tNew.Alleles\n")

if __name__ == "__main__":
    main()