#!/usr/bin/env python
import sys
from collections import defaultdict
import pandas as pd


def main():
    infile, outfile = sys.argv[1:]
    d = pd.read_csv(infile, sep="\t" if infile.endswith(".tsv") else ",", index_col=0)
    d = d[d["Size"] >= 2]
    counter = defaultdict(int)
    for t, tc in d[["T", "T-C"]].values:
        if tc < 2:
            tc = 0
        counter[(t, tc)] += 1
    
    m = pd.DataFrame([[t, tc, count] for (t, tc), count in sorted(counter.items()) ])
    m.columns = ["#T", "TC", "Count"]
    m.to_csv(outfile, sep="\t" if outfile.endswith(".tsv") else ",", index=False)


if __name__ == "__main__":
    main()