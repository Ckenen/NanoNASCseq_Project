#!/usr/bin/env python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd

def main():

    infile, outfile = sys.argv[1:]

    d = pd.read_csv(infile, sep="\t", index_col=0)
    d = d[d["Size"] >= 2]
    base_counter = defaultdict(int)
    for base in "ACGT":
        base_counter[base] = d[base].sum()
    event_counter = defaultdict(int)
    for ref in "ACGT":
        for alt in "ACGT":
            if ref == alt:
                continue
            s = d["%s-%s" % (ref, alt)]
            event_counter["%s%s" % (ref, alt)] = sum(s[s >= 2])

    ks = []
    rows = []
    for k, v in sorted(event_counter.items()):
        ks.append(k)
        rows.append([k[0], k[1], base_counter[k[0]], v, np.divide(v, base_counter[ref])])
    m = pd.DataFrame(rows, index=pd.Index(ks, name="Type"), columns=["RefBase", "AltBase", "BaseCount", "EventCount", "Ratio"])
    m.to_csv(outfile, sep="\t" if outfile.endswith(".tsv") else ",")
    
    
if __name__ == "__main__":
    main()