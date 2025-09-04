#!/usr/bin/env python
import sys
import pandas as pd
import anndata as ad

# single-cell count matrix

def main():
    f_filelist, outfile = sys.argv[1:]
    
    array1 = []
    array2 = []
    for path in open(f_filelist):
        path = path.strip()
        assert path[-4:] == ".tsv"
        cell = path.split("/")[-1][:-4]
        d = pd.read_csv(path, sep="\t", header=0, index_col=0)
        s1 = d["Total"]
        s2 = d["New"]
        s1.name = cell
        s2.name = cell
        array1.append(s1)
        array2.append(s2)    

    dat1 = pd.DataFrame(array1).fillna(0)
    dat2 = pd.DataFrame(array2).fillna(0)
    dat1.index.name = "Cell"
    dat2.index.name = "Cell"
    dat1 = dat1[sorted(dat1.columns)] # sort gene id
    dat2 = dat2[sorted(dat2.columns)]
    assert all(dat1.index == dat2.index)
    assert all(dat1.columns == dat2.columns)

    adata = ad.AnnData(dat1)
    adata.layers["total"] = dat1
    adata.layers["new"] = dat2
    adata.write(outfile, compression="gzip")
    

if __name__ == "__main__":
    main()
