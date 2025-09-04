#!/usr/bin/env python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import pysam


def main():
    vcffile, bamfile, outfile = sys.argv[1:]
    
    rows = []
    with pysam.VariantFile(vcffile) as vcf, pysam.AlignmentFile(bamfile) as bam:
        for record in vcf:
            segments = defaultdict(list)
            for s in bam.fetch(record.chrom, record.start, record.stop):
                segments[s.query_name].append(s)
            hp0 = 0
            hp1 = 0
            hp2 = 0
            for name, ss in segments.items():
                if len(ss) == 1:
                    s = ss[0]
                    if s.has_tag("HP"):
                        hp = s.get_tag("HP")
                        if hp == 1:
                            hp1 += 1
                        elif hp == 2:
                            hp2 += 1
                        else:
                            assert False
                    else:
                        hp0 += 1
                elif len(ss) == 2:
                    s1, s2 = ss
                    if s1.has_tag("HP") and s2.has_tag("HP"):
                        v1, v2 = s1.get_tag("HP"), s2.get_tag("HP")
                        if v1 == v2:
                            if v1 == 1:
                                hp1 += 1
                            elif v1 == 2:
                                hp2 += 1
                            else:
                                assert False
                        else:
                            hp0 += 1
                    else:
                        hp0 += 1
                else:
                    assert False
            row = [record.chrom, record.start, record.stop, hp0, hp1, hp2]
            rows.append(row)
            
    dat = pd.DataFrame(rows)
    dat.columns = ["Chrom", "Start", "End", "UnAssinged", "HP1", "HP2"]
    dat["Assigned"] = dat["HP1"] + dat["HP2"]
    dat["Assigned.Ratio"] = dat["Assigned"] / (dat[["UnAssinged", "HP1", "HP2"]].sum(axis=1))
    dat["Log2FC"] = np.log2(dat["HP1"] / dat["HP2"])
    dat.to_csv(outfile, sep="\t", index=False)
    

if __name__ == "__main__":
    main()