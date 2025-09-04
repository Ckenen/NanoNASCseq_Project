#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt


def main():
    infile, outfile = sys.argv[1:]
    
    name = os.path.basename(outfile)
    dat = pd.read_csv(infile, sep="\t", index_col=0)
    dat["Ratio"] = dat["Reads"] / dat["Reads"].sum()    
    cells = list(filter(lambda item: item != "unclassified", dat.index))
    xticks = [c.split(".")[-1] for c in cells]
    tmp = dat.loc[cells]
    ys = tmp["Reads"] / 1e6
    xs = np.arange(len(ys))
    
    plt.figure(figsize=((len(cells) + 6) * 0.25, 3))
    plt.title(name)
    plt.bar(xs, ys, color="C0", edgecolor="black")
    plt.xlim(min(xs) - 0.5, max(xs) + 0.5)
    plt.xticks(xs, xticks, rotation=90)
    plt.ylabel("Read number (M)")
    plt.grid(axis="y", ls="--")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()
    
    
if __name__ == '__main__':
    main()
    