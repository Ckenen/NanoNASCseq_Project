#!/usr/bin/env python
import sys
import pandas as pd


def main():
    infile = sys.argv[1]
    categories = ['full-splice_match', 
                  'incomplete-splice_match', 
                  'novel_in_catalog', 
                  'novel_not_in_catalog']
    dat = pd.read_csv(infile, sep="\t", index_col=0)
    dat = dat[[x in categories for x in dat["structural_category"]]]
    genes = len(set(dat["associated_gene"]))
    cell = infile.split("/")[-2]
    print(cell, genes, sep="\t")
            
    
if __name__ == "__main__":
    main()