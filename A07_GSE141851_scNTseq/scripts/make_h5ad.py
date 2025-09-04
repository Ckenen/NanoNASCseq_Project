#!/usr/bin/env python
import sys
import glob
import os
from collections import defaultdict
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad 

def main():    
    infile, outfile = sys.argv[1:]
    
    data_total = defaultdict(int)
    data_nascent = defaultdict(int)
    cell_umis = defaultdict(int)
    gene_umis = defaultdict(int)
    with open(infile) as f:
        for line in f:
            name, cell, umi, reads = line.strip("\n").split("\t")
            suffix = name[-3:]
            gene = name[:-3]
            if suffix == "--T" or suffix == "--C":
                data_total[(cell, gene)] += 1
                if suffix == "--C":
                    data_nascent[(cell, gene)] += 1
                cell_umis[cell] += 1
                gene_umis[gene] += 1
                
    cells_filtered = set()
    for cell, umis in cell_umis.items():
        if umis >= 500:
            cells_filtered.add(cell)
    # print(len(cells_filtered))

    genes_filtered = set()
    for gene, umis in gene_umis.items():
        if umis >= 10:
            genes_filtered.add(gene)
    # print(len(genes_filtered))

    genes_filtered2 = set()
    for (cell, gene), umis in data_total.items():
        if cell in cells_filtered and gene in genes_filtered:
            genes_filtered2.add(gene)
    # print(len(genes_filtered2))

    cells = list(sorted(cells_filtered))
    genes = list(sorted(genes_filtered))
    rows1 = []
    rows2 = []
    for cell in cells:
        rows1.append([data_total[(cell, gene)] for gene in genes])
        rows2.append([data_nascent[(cell, gene)] for gene in genes])

    adata = ad.AnnData(pd.DataFrame(rows1, index=pd.Index(cells, name="Cell"), columns=genes))
    adata.layers["total"] = np.matrix(rows1)
    adata.layers["nascent"] = np.matrix(rows2)
    adata.var.index = pd.Index(genes, name="Gene")
    adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
    sc.pp.filter_cells(adata, min_genes=500)
    sc.pp.filter_cells(adata, max_genes=5000)
    sc.pp.filter_genes(adata, min_cells=10)
    adata.write_h5ad(outfile, compression="gzip")

if __name__ == "__main__":
    main()