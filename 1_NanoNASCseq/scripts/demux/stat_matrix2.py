#!/usr/bin/env python
import sys
import gzip
from collections import defaultdict


def main():
    infile, outfile = sys.argv[1:]
    
    data = defaultdict(int)
    
    if infile.endswith(".gz"):
        f = gzip.open(infile, "rt")
    else:
        f = open(infile)
    for line in f:
        row = line.strip("\n").split("\t")
        if int(row[1]) < 400:
            continue
        data[(row[2], row[3], row[4], row[7], row[8], row[9], row[10], row[13])] += 1
    f.close()
    
    if outfile.endswith(".gz"):
        fw = gzip.open(outfile, "wt")
    else:
        fw = open(outfile, "w+")
    fw.write("Barcode1,Direction1,Location1,ED1,Barcode2,Direction2,Location2,ED2,Count\n")
    for k, c in sorted(data.items()):
        row = list(k)
        row.append(c)
        fw.write(",".join(map(str, row)) + "\n")
    fw.close()
    
    
if __name__ == "__main__":
    main()