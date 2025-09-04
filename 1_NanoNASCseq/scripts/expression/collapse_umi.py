#!/usr/bin/env python
import sys
from collections import defaultdict, Counter
import numpy as np
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import BamFile
from pyBioInfo.Utils import BundleBuilder, BlockTools


def collapse_edge(positions, edge_type):
    counts = Counter(positions).items() # [position, count]
    counts = list(sorted(counts, key=lambda x: (x[1], x[0]), reverse=True))

    # unique position or maximum support position
    if len(counts) == 1 or counts[0][1] > counts[1][1]:
        return counts[0][0]
    
    # closest to mean position
    max_count = counts[0][1]
    counts1 = list(filter(lambda x: x[1] == max_count, counts)) # [position, max_count]
    mean = np.mean([x[0] for x in counts1])
    distances = [[x[0], abs(x[0] - mean)] for x in counts1]
    distances = list(sorted(distances, key=lambda x: x[1]))
    min_distance = distances[0][1]
    distances1 = list(filter(lambda x: x[1] == min_distance, distances))
    if len(distances1) == 1:
        return distances1[0][0]

    # for start, select minimum position
    # for end, select maximum position
    ps = [x[0] for x in distances1]
    ps.sort()
    if edge_type == "start":
        return ps[0]
    elif edge_type == "end":
        return ps[-1]
    assert False


def collapse_introns(cluster_name, reads, start, end):
    # assert len(reads) > 0
    # if len(reads) == 1:
    #     pass
    # elif len(reads) == 2:
    #     pass
    # else:
    #     pass

    for read in reads:
        read.introns = set([(x, y) for x, y in BlockTools.gaps(read.blocks)])
        
    introns = []
    for read in reads:
        for intron in read.introns:
            introns.append(intron)
    array = list(sorted(Counter(introns).items()))
    # print(array)
    
    i = 0
    while i < len(array) - 1:
        intron1, count1 = array[i]
        intron2, count2 = array[i + 1]
        x = max(intron1[0], intron2[0])
        y = min(intron1[1], intron2[1])
        if x <= y:
            if count1 >= count2:
                array.pop(i + 1)
            elif count1 < count2:
                array.pop(i)
            else:
                sys.stderr.write("%s\n" % reads[0].format("bed"))
                sys.stderr.write("%s, %s\n" % (intron1, intron2))
                raise RuntimeError()
        else:
            i += 1
    
    tmp = []
    for intron in array:
        x, y = intron[0]
        if x > start and y < end:
           tmp.append(intron) 
        else:
            continue
            # print(intron, start, end, sep="\t")
            # assert False
    array = tmp
    
    finals = []
    for intron in array:
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
        ratio = agree / covered
        if ratio >= 0.5:
            finals.append(pos)
        elif ratio < 0.5:
            pass
        else:
            print(cluster_name)
            print(pos, support, covered, agree, disagree, ratio, sep="\t")
            assert False
    return finals


def main():
    default_cn = 0
    with BamFile(sys.argv[1], random=True) as f:
        for bundle in BundleBuilder(f, keep=True):
            clusters = defaultdict(list)
            for read in bundle.data:
                if read.segment.has_tag("CN"):
                    cn = read.segment.get_tag("CN") # cluster name
                else:
                    cn = "UMI.%d" % default_cn
                    default_cn += 1
                clusters[cn].append(read)
            for cluster_name, reads in clusters.items():
                try:                          
                    start = collapse_edge([read.start for read in reads], edge_type="start")
                    end = collapse_edge([read.end for read in reads], edge_type="end")
                    introns = collapse_introns(cluster_name, reads, start, end)
                    blocks = []
                    p1 = start
                    for intron in introns:
                        assert p1 < intron[0]
                        blocks.append([p1, intron[0]])
                        p1 = intron[1]
                    assert p1 < end
                    blocks.append([p1, end])
                    obj = GRange(chrom=reads[0].chrom, 
                                blocks=blocks, 
                                name=cluster_name, 
                                strand=reads[0].strand)
                    obj.score = len(reads)
                    # obj.score = reads[0].segment.get_tag("CS")
                    print(obj.format("BED"))
                except Exception:
                    continue

                
if __name__ == '__main__':
    main()
    
