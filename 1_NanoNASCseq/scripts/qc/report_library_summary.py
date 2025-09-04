#!/usr/bin/env python
import sys
import os
import gzip
import pandas as pd


def main():
    infile = sys.argv[1]
    
    run = os.path.splitext(os.path.basename(infile))[0]

    dat = pd.read_csv(infile, sep="\t", header=0)
    dat = dat.sort_values(by="Length")


    bases = 0
    reads = 0
    for length, count in dat[["Length", "Count"]].values:
        bases += length * count
        reads += count
    mean = bases / reads

    median = None
    if reads % 2 == 1:
        i = int(reads / 2)
        j = 0
        for length, count in dat[["Length", "Count"]].values:
            j += count
            if j > i:
                median = length
                break
    else:
        i2 = int(reads / 2)
        i1 = i2 - 1
        j = 0
        v1 = None
        for length, count in dat[["Length", "Count"]].values:
            j += count
            if j > i1 and v1 is None:
                v1 = length
            if j > i2:
                median = (length + v1) / 2
                break
                
    n50 = None
    i = bases / 2
    j = 0
    for length, count in dat[["Length", "Count"]].values[::-1]:
        j += length * count
        if j >= i:
            n50 = length
            break
        
    print("Run", "TotalBase", "TotalRead", "MaxLength", "MeanLength", "MedianLength", "N50", sep="\t")
    print(run, bases, reads, max(dat["Length"]), mean, median, n50, sep="\t")


if __name__ == '__main__':
    main()
