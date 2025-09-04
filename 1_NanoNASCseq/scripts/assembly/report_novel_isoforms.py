#!/usr/bin/env python
import sys
import pandas as pd
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder


def main():
    f_tsv, f_gtf = sys.argv[1:]

    d = pd.read_csv(f_tsv, sep="\t", header=0, index_col=0)
    sc_list = ["novel_in_catalog", "novel_not_in_catalog"]
    d = d[[x in sc_list for x in d["structural_category"]]]
    isoform_list = set(d.index)

    with GtfFile(f_gtf) as f:
        records = [x for x in f]
        isoforms = list(GtfTranscriptBuilder(records))

    tmp = []
    for x in isoforms:
        if x.name in isoform_list:
            tmp.append(x)
    isoforms = tmp

    for x in isoforms:
        for k, v in x.records.items():
            for r in v:
                print(r.format())


if __name__ == "__main__":
    main()