#!/usr/bin/env python
import sys
import os
import glob
from collections import defaultdict
import pandas as pd
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder, BedFile
from pyBioInfo.Utils import BlockTools


def cluster_by_cap(clusters):
    new_clusters = []
    for ts in clusters:
        ts = list(sorted(ts, key=lambda t: t.cap))
        cluster = None
        for t in ts:
            if cluster is None:
                cluster = [t]
            elif t.cap - cluster[-1].cap <= 500:
                cluster.append(t)
            else:
                new_clusters.append(cluster)
                cluster = [t]
        if cluster is not None:
            new_clusters.append(cluster)
    return new_clusters


def cluster_by_polya(clusters):
    new_clusters = []
    for ts in clusters:
        ts = list(sorted(ts, key=lambda t: t.polya))
        cluster = None
        for t in ts:
            if cluster is None:
                cluster = [t]
            elif t.polya - cluster[-1].polya <= 200:
                cluster.append(t)
            else:
                new_clusters.append(cluster)
                cluster = [t]
        if cluster is not None:
            new_clusters.append(cluster)
    return new_clusters


def main():
    f_config, outdir = sys.argv[1:]

    config = pd.read_csv(f_config, sep="\t", header=0)

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if not os.path.exists(outdir + "/source"):
        os.mkdir(outdir + "/source")

    fws = dict()
    for cell, f_gtf in config[["Cell", "Gtf"]].values:
        with GtfFile(f_gtf) as f:
            records = [x for x in f]
            transcripts = list(sorted(GtfTranscriptBuilder(records)))
            for t in transcripts:
                chrom = t.chrom
                if chrom not in fws:
                    outfile = outdir + "/source/%s.bed" % chrom
                    fws[chrom] = open(outfile, "w+")
                fw = fws[chrom]
                t.name = "%s:%s" % (cell, t.name)
                fw.write(t.format("bed") + "\n")
    for fw in fws.values():
        fw.close()

    array = []
    for f_bed in sorted(glob.glob(outdir + "/source/*.bed")):
        transcripts = defaultdict(list)
        with BedFile(f_bed) as f:
            for t in f:
                if len(t.blocks) > 1:
                    t.introns = tuple(BlockTools.gaps(t.blocks))
                    t.cell = t.name.split(":")[0]
                    if t.strand == "+":
                        t.cap = t.start
                        t.polya = t.end
                    else:
                        t.cap = t.end
                        t.polya = t.start
                    k = (t.chrom, t.introns, t.strand)
                    transcripts[k].append(t)

        for k, v in transcripts.items():
            chrom, introns, strand = k
            clusters = [v]
            repeat = 0
            num_cluster = 0
            while True:
                clusters = cluster_by_polya(clusters)
                clusters = cluster_by_cap(clusters)
                if repeat == 0 or len(clusters) != num_cluster:
                    num_cluster = len(clusters)
                else:
                    break
                repeat += 1

            for cluster in clusters:
                support_cells = len(set([t.cell for t in cluster]))
                start = min([t.start for t in cluster])
                end = max([t.end for t in cluster])
                blocks = []
                s = start
                for intron in introns:
                    blocks.append([s, intron[0]])
                    s = intron[1]
                blocks.append([s, end])
                obj = GRange(chrom=chrom, blocks=blocks, name="Unknown", strand=strand)
                obj.score = support_cells
                array.append(obj)

    # minimum support by 5 cells
    array = list(filter(lambda t: t.score >= 5, array))
    array.sort()
    
    for i, t in enumerate(array):
        t.name = "Merged.%d" % i

    with open(outdir + "/merged.bed", "w+") as fw:
        for t in array:
            fw.write(t.format("bed") + "\n")

    with open(outdir + "/merged.gtf", "w+") as fw:
        for t in array:
            s = 'gene_id "%s"; transcript_id "%s"; support_cells %d;' % (t.name, t.name, t.score)
            # transcript
            row = [t.chrom, "Merged", "transcript", t.start + 1, t.end, ".", t.strand, ".", s]
            line = "\t".join(map(str, row))
            fw.write(line + "\n")
            # exon
            for exon in t.blocks:
                row = [t.chrom, "Merged", "exon", exon[0] + 1, exon[1], ".", t.strand, ".", s]
                line = "\t".join(map(str, row))
                fw.write(line + "\n")


if __name__ == "__main__":
    main()