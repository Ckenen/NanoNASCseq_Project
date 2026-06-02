#!/usr/bin/env python
import sys, os, glob, pickle
import pandas as pd
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
import seaborn as sns


def main():
    # SC_EXP_FNAME = "GSE60361_C1-3005-Expression.txt"
    # OUTDIR = "./"
    # SC_EXP_FNAME = "results/blastocyst/mouse_gene_counts.tsv"
    # GENOME = "mm10"
    # OUTDIR = "results/blastocyst/pySCENIC_output"
    
    SC_EXP_FNAME, GENOME, OUTDIR = sys.argv[1:]
    
    if not os.path.exists(OUTDIR):
        os.mkdir(OUTDIR)
        
    SCENIC_ROOT = "~/software/SCENIC" # TODO: change to your SCENIC cisTarget databases root directory
    
    # allTFs_*.txt are downloaded from https://resources.aertslab.org/cistarget/tf_lists/
    # *.feather are downloaded from https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/
    # and https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/
    # motifs*.tbl are downloaded from https://resources.aertslab.org/cistarget/motif2tf/
    # You should download the necessary files from these locations to the SCENIC_ROOT directory!
    
    if GENOME == "mm10":
        MM_TFS_FNAME = os.path.join(SCENIC_ROOT, "allTFs_mm.txt")
        DATABASES_GLOB = os.path.join(SCENIC_ROOT, "mm10_*.feather")
        MOTIF_ANNOTATIONS_FNAME = os.path.join(SCENIC_ROOT, "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
    elif GENOME == "hg38":
        MM_TFS_FNAME = os.path.join(SCENIC_ROOT, "allTFs_hg38.txt")
        DATABASES_GLOB = os.path.join(SCENIC_ROOT, "hg38_*.feather")
        MOTIF_ANNOTATIONS_FNAME = os.path.join(SCENIC_ROOT, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
    else:
        raise NotImplementedError("Genome not supported: {}".format(GENOME))
        # assert False
    
    ADJACENCIES_FNAME = os.path.join(OUTDIR, "adjacencies.csv")
    MODULES_FNAME = os.path.join(OUTDIR, "modules.p")
    REGULONS_FNAME = os.path.join(OUTDIR, "regulons.p")
    MOTIFS_FNAME = os.path.join(OUTDIR, "motifs.csv")
    AUC_FNAME = os.path.join(OUTDIR, "auc.csv")
    BINARY_FNAME = os.path.join(OUTDIR, "binary.csv")
    HEATMAP_FNAME = os.path.join(OUTDIR, "auc_heatmap.pdf")

    # input table: rowname is gene, colname is cell
    # change to: rowname is cell, colname is gene
    ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
    
    tf_names = load_tf_names(MM_TFS_FNAME)

    db_fnames = glob.glob(DATABASES_GLOB)
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

    # Adjacencies
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
    adjacencies.to_csv(ADJACENCIES_FNAME, sep="\t", index=False)
    
    # Modules
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    with open(MODULES_FNAME, 'wb') as f:
        pickle.dump(modules, f)
        
    # Calculate a list of enriched motifs and the corresponding target genes for all modules.
    with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)
    # Save the enriched motifs and the discovered regulons to disk.
    df.to_csv(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "wb") as f:
        pickle.dump(regulons, f)
    
    # AUC
    auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
    auc_mtx.to_csv(AUC_FNAME, index=True, sep='\t')
    ret = sns.clustermap(auc_mtx, figsize=(8,8))
    ret.fig.savefig(HEATMAP_FNAME, dpi=300)
    
    # Binary
    binary_mtx = binarize(auc_mtx)
    binary_mtx[0].to_csv(BINARY_FNAME, index=True, sep='\t')
    # plot_binarization(auc_mtx,auc_mtx.columns[0])
    
    
if __name__ == "__main__":
    main()