#!/usr/bin/env python
import sys
import pandas as pd

# single-cell count matrix

def main():
    f_filelist, out_total, out_nascent = sys.argv[1:]
    
    array1 = []
    array2 = []

    for path in open(f_filelist):
        path = path.strip()
        name = path.split("/")[-1][:-4]
        d = pd.read_csv(path, sep="\t", header=0, index_col=0)
        index_name = d.index.name
        
        s1 = d["Total"]
        s1.name = name
        array1.append(s1)
        
        s2 = d["Nascent"]
        s2.name = name
        array2.append(s2)
        
    dat1 = pd.concat(array1, axis=1).fillna(0)
    dat1.index.name = index_name
    dat1 = dat1.sort_index()
    dat2 = pd.concat(array2, axis=1).fillna(0)
    dat2.index.name = index_name
    dat2 = dat2.sort_index()
    
    dat1.to_csv(out_total, sep="\t")
    dat2.to_csv(out_nascent, sep="\t")
    

if __name__ == "__main__":
    main()