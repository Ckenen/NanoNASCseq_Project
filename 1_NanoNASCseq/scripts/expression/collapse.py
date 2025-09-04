#!/usr/bin/env python
import sys
import os
import optparse
from collections import defaultdict, Counter
import numpy as np
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import BamFile
from pyBioInfo.Utils import BundleBuilder, BlockTools


def determine_edge(positions, select_max=True):
    counts = Counter(positions)
    items = list(sorted(counts.items(), key=lambda item: item[1], reverse=True))

    if len(items) == 1 or items[0][1] > items[1][1]:
        return items[0][0]
    
    ps = []
    for item in items:
        if item[1] == items[0][1]:
            ps.append(item[0])
    assert len(ps) > 1
    
    mean = np.mean(positions)
    distances = []
    for p in ps:
        d = abs(p - mean)
        distances.append([p, d])
    distances = list(sorted(distances, key=lambda item: item[1]))
    
    ps = []
    for p, d in distances:
        if d == distances[0][1]:
            ps.append(p)
    assert len(ps) >= 1
    if len(ps) == 1:
        return ps[0]
    ps.sort()
    if select_max:
        return ps[-1]
    else:
        return ps[0]
    # assert False      


def determine_introns(reads, start, end):
    for read in reads:
        read.introns = set([(x, y) for x, y in BlockTools.gaps(read.blocks)])
        
    introns = []
    for read in reads:
        for intron in read.introns:
            introns.append(intron)
    counts = Counter(introns)
    introns = sorted(counts.items())
    
    i = 0
    while i < len(introns) - 1:
        intron1, intron2 = introns[i], introns[i + 1]
        x = max(intron1[0][0], intron2[0][0])
        y = min(intron1[0][1], intron2[0][1])
        if x <= y:
            c1, c2 = intron1[1], intron2[1]
            if c1 >= c2:
                introns.pop(i + 1)
            elif c1 < c2:
                introns.pop(i)
            else:
                print(intron1, intron2)
                assert False
        else:
            i += 1
    
    tmp = []
    for intron in introns:
        x, y = intron[0]
        if x > start and y < end:
           tmp.append(intron) 
        else:
            continue
            # print(intron, start, end, sep="\t")
            # assert False
    introns = tmp
    
    finals = []
    for intron in introns:
        pos = intron[0]
        support = intron[1]
        covered = 0
        agree = 0
        disagree = 0
        
        for read in reads:
            if read.start < pos[1] and read.end > pos[0]:
                covered += 1
                if pos in read.introns:
                    agree += 1
                else:
                    disagree += 1
        # assert agree == support
        ratio = agree / covered
        if ratio >= 0.5:
            finals.append(pos)
        elif ratio < 0.5:
            pass
        else:
            print(reads[0].segment.get_tag("CN"))
            print(pos, support, covered, agree, disagree, ratio, sep="\t")
            assert False
        # print(pos, support, covered, agree, disagree, sep="\t")
    return finals


def main():
    parser = optparse.OptionParser(usage="%prog [options] input.bam")
    parser.add_option("-b", "--bed", dest="bed", metavar="PATH", help="[%default]")
    parser.add_option("-g", "--gtf", dest="gtf", metavar="PATH", help="[%default]")
    options, args = parser.parse_args()
    f_bam = args[0]
    f_bed = options.bed
    f_gtf = options.gtf
    
    # bamfile, prefix = sys.argv[1:]
    # f_gtf = prefix + ".unsorted.gtf"
    # f_sorted_gtf = prefix + ".gtf.gz"
    # f_bed = prefix + ".unsorted.bed"
    # f_sorted_bed = prefix + ".bed.gz"
        
    fw_gtf = open(f_gtf, "w+") if f_gtf else None
    fw_bed = open(f_bed, "w+") if f_bed else None
        
    with BamFile(f_bam) as f:
        for bundle in BundleBuilder(f, keep=True):
            clusters = defaultdict(list)
            for read in bundle.data:
                cluster_name = read.segment.get_tag("CN") # umi cluster name
                clusters[cluster_name].append(read)
                
            for cluster_name, reads in clusters.items():                                
                chrom = reads[0].chrom
                strand = reads[0].strand
                start = determine_edge([read.start for read in reads], select_max=False)
                end = determine_edge([read.end for read in reads], select_max=True)
                # print(cluster_name, chrom, start, end, sep="\t")
        
                introns = determine_introns(reads, start, end)
                # print(introns)
                exons = []
                p1, p2 = start, end
                for intron in introns:
                    exons.append([p1, intron[0]])
                    p1 = intron[1]
                exons.append([p1, p2])
                # print(exons)
                
                if fw_gtf:
                    gid = "GID.%s" % cluster_name
                    tid = "TID.%s" % cluster_name
                    attris = "gene_id \"%s\"; transcript_id \"%s\"; umi_size \"%d\";" % (gid, tid, len(reads))
                    line = "\t".join(map(str, [chrom, "RNA", "transcript", 
                                        start + 1, end, ".", 
                                        strand, ".", attris]))
                    fw_gtf.write(line + "\n")
                    for exon in exons:
                        line = "\t".join(map(str, [chrom, "RNA", "exon", 
                                            exon[0] + 1, exon[1], ".", 
                                            strand, ".", attris]))
                        fw_gtf.write(line + "\n")
                    
                if fw_bed:
                    obj = GRange(chrom=chrom, start=start, end=end, name=cluster_name, strand=strand, blocks=exons)
                    fw_bed.write(obj.format("BED") + "\n")
                    
    if fw_gtf:
        fw_gtf.close()
    if fw_bed:
        fw_bed.close()
  
                
if __name__ == '__main__':
    main()
    
