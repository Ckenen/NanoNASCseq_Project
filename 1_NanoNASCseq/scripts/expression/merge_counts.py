#!/usr/bin/env python
import sys
from collections import defaultdict
import pandas as pd

# pseudo-bulk expression
            
def main():
    filelist, outfile = sys.argv[1:]

    paths = [line.strip() for line in open(filelist)]
    
    data = dict()
    header = None
    for path in paths:
        m = pd.read_csv(path, sep="\t", header=0)
        header = list(m.columns)
        for gid, total, new, total_alleles, new_alleles in m.values:
            if gid not in data:
                data[gid] = {
                    "Total": 0, 
                    "New": 0, 
                    "Total.Alleles": defaultdict(int), 
                    "New.Alleles": defaultdict(int)}
            data[gid]["Total"] += total
            data[gid]["New"] += new
            if isinstance(total_alleles, str):
                for item in total_alleles.split(";"):
                    allele, count = item.split(":")
                    data[gid]["Total.Alleles"][allele] += int(count)
            if isinstance(new_alleles, str):
                for item in new_alleles.split(";"):
                    allele, count = item.split(":")
                    data[gid]["New.Alleles"][allele] += int(count)
    header.append("Cells")
    
    sep = ","
    if outfile.endswith(".tsv"):
        sep = "\t"
    with open(outfile, "w+") as fw:
        fw.write(sep.join(header) + "\n")
        for gid in sorted(data.keys()):
            total = data[gid]["Total"]
            new = data[gid]["New"]
            total_alleles = data[gid]["Total.Alleles"]
            new_alleles = data[gid]["New.Alleles"]
            s1 = ";".join(["%s:%d" % (k, v) for k, v in sorted(total_alleles.items())])
            s2 = ";".join(["%s:%d" % (k, v) for k, v in sorted(new_alleles.items())])
            line = sep.join(map(str, [gid, total, new, s1, s2, len(paths)]))
            fw.write(line + "\n")
        

if __name__ == "__main__":
    main()
