#!/usr/bin/env python
import os, optparse
from collections import defaultdict
import pandas as pd


def make_one_row_mismatch_ratio(infile, outfile):
    d = pd.read_csv(infile, sep="\t" if infile.endswith(".tsv") else ",", index_col=0)
    bases = ["A", "C", "G", "T"]
    base_counter = defaultdict(int)
    for base in bases:
        count = d[d["RefBase"] == base]["BaseCount"].values[0]
        base_counter[base] = count
    event_counter = defaultdict(int)
    event_ratios = defaultdict(float)
    for ref in ["A", "C", "G", "T"]:
        for alt in ["A", "C", "G", "T"]:
            if ref == alt:
                continue
            e = "%s%s" % (ref, alt)
            if "EventCount.NoSNP" in d.columns:
                count = d["EventCount.NoSNP"][e]
            else:
                count = d["EventCount"][e]
            event_counter[e] = count
            event_ratios[e + "_ratio"] = count / base_counter[ref]
    name = os.path.splitext(os.path.basename(infile))[0]
    header = ["Name"] + list(base_counter.keys()) + list(event_counter.keys()) + list(event_ratios.keys())
    row = [name] + list(base_counter.values()) + list(event_counter.values()) + list(event_ratios.values())
    m = pd.DataFrame([row], columns=header)
    m.to_csv(outfile, sep="\t" if outfile.endswith(".tsv") else ",", index=False)


if __name__ == "__main__":
    parser = optparse.OptionParser(usage="%prog -i ratios.csv -o ratios.csv")
    parser.add_option("-i", "--infile", dest="infile")
    parser.add_option("-o", "--outfile", dest="outfile")
    options, _ = parser.parse_args()
    make_one_row_mismatch_ratio(
        infile=options.infile,
        outfile=options.outfile,
    )
