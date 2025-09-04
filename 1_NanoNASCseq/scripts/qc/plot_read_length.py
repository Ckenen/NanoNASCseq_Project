#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    infile, outfile = sys.argv[1:]

    name = os.path.splitext(os.path.basename(outfile))[0]
    dat = pd.read_csv(infile, sep="\t", header=0)
    dat = dat.sort_values(by="Length")
    
    reads = sum(dat["Count"])
    bases = sum(dat["Length"] * dat["Count"])
    
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
            
    MIN_LENGTH = 0
    MAX_LENGTH = 5000
    BIN_WIDTH = 25
    BIN_COUNT = int((MAX_LENGTH - MIN_LENGTH) / BIN_WIDTH)
    ys = np.zeros(BIN_COUNT, dtype=int)
    xs = np.arange(len(ys)) * BIN_WIDTH + BIN_WIDTH / 2
    for length, count in dat[["Length", "Count"]].values:
        i = int(length / BIN_WIDTH)
        if 0 <= i < len(ys):
            ys[i] += count
    ys = ys / BIN_WIDTH

    plt.figure(figsize=(5, 3))
    plt.title(name)
    plt.plot(xs, ys, color="C0", lw=1)
    plt.axvline(n50, lw=1, ls="--", color="C1")
    plt.axvline(mean, lw=1, ls="--", color="C2")
    plt.axvline(median, lw=1, ls="--", color="C3")
    # plt.axhline(0, lw=1, color="grey")
    max_y = max(ys)
    plt.text(MAX_LENGTH * 0.5, max_y * 0.95, "Reads=%s" % format(reads, ","))
    plt.text(MAX_LENGTH * 0.5, max_y * 0.85, "Bases=%s" % format(bases, ","))
    plt.text(MAX_LENGTH * 0.5, max_y * 0.75, "N50=%s" % format(round(n50, 2), ","), color="C1")
    plt.text(MAX_LENGTH * 0.5, max_y * 0.65, "Mean=%s" % format(round(mean, 2), ","), color="C2")
    plt.text(MAX_LENGTH * 0.5, max_y * 0.55, "Median=%s" % format(round(median, 2), ","), color="C3")
    plt.xlabel("Read length (nt)")
    plt.ylabel("Read numbers")
    plt.grid(ls="--", lw=1, color="lightgrey")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
