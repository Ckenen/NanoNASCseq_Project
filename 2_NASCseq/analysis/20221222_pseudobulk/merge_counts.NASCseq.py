#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd


def get_halflife(total, nascent, time):
    if time == 0 or total == 0:
        return np.nan
    elif nascent == 0:
        return np.inf
    elif total == nascent:
        return 0
    else:
        return -time/np.log2(1-nascent/total)


def main():
    filelist, annofile, time, outfile = sys.argv[1:]
    time = float(time) # hours
    
    paths = [line.strip() for line in open(filelist)]
    
    m = None
    for path in paths:
        d = pd.read_csv(path, sep="\t", index_col=0)
        d = d[["Length", "Count", "Count.Old", "Count.New", "LibSize"]]
        if m is None:
            m = d.copy()
        else:
            for c in ["Count", "Count.Old", "Count.New", "LibSize"]:
                m[c] = m[c] + d[c]
    m["FPKM"] = m["Count"] * 1e9 / m["Length"] / m["LibSize"]
    m["FPKM.Old"] = m["Count.Old"] * 1e9 / m["Length"] / m["LibSize"]
    m["FPKM.New"] = m["Count.New"] * 1e9 / m["Length"] / m["LibSize"]
    m["Cells"] = len(paths)
    m["NTR"] = m["Count.New"] / m["Count"]
    m["Halflife"] = [get_halflife(total, nascent, time) for total, nascent in m[["Count", "Count.New"]].values]    

    anno = pd.read_csv(annofile, index_col=0)
    anno = anno.drop(columns=["Length"])
    m = m.merge(anno, left_index=True, right_index=True)
    m.to_csv(outfile, sep="\t" if outfile.endswith(".tsv") or outfile.endswith(".tsv.gz") else ",")


if __name__ == "__main__":
    main()