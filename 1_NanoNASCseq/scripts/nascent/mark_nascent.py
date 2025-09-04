#!/usr/bin/env python
import sys
from collections import defaultdict
import numpy as np
from pyBioInfo.IO.File import Alignment
from pyBioInfo.Utils import BundleBuilder
import pysam


def load_alignments(f, chrom):
    for segment in f.fetch(chrom):
        obj = Alignment(segment)
        
        events = []
        for item in obj.segment.get_tag("CE").split(";"):
            if item == "":
                continue
            e = item.split(",")
            e[0], e[4] = int(e[0]), int(e[4])
            if "/" in e[3]:
                e[3] = [int(x) for x in e[3].split("/")]
            else:
                e[3] = int(e[3])
            events.append(e)
        obj.ce = events
        
        refs = defaultdict(int)
        for item in obj.segment.get_tag("RC").split(";"):
            ref, count = item.split(",")
            count = int(count)
            refs[ref] = count
            
        tc_locs = []
        if obj.strand == "+":
            ref, alt = "T", "C"
        else:
            ref, alt = "A", "G"
        for e in obj.ce:
            if e[1] == ref and e[2] == alt:
                tc_locs.append(e[0])
                
        obj.tc_locs = tc_locs
        obj.tc_count = len(tc_locs)
        obj.t_count = refs[ref]
        obj.tc_ratio = np.divide(obj.tc_count, obj.t_count)
        obj.cluster_name = segment.get_tag("CN") # umi cluster name
        
        yield obj


def main():
    infile, outfile, outfile2 = sys.argv[1:]
    
    f = pysam.AlignmentFile(infile)
    fw = pysam.AlignmentFile(outfile, "wb", f)
    fw2 = open(outfile2, "w+")
    
    for chrom in f.references:
        print("Processing %s" % chrom)
        
        for bundle in BundleBuilder(load_alignments(f, chrom), keep=True):
            clusters = defaultdict(list)
            for obj in bundle.data:
                clusters[obj.cluster_name].append(obj)
                
            for name, ds in clusters.items():
                # ps = get_ds_signals(ds)
                counter = defaultdict(int)
                for obj in ds:
                    for p in obj.tc_locs:
                        counter[p] += 1
                s = ";".join(["%d:%d" % (k, v) for k, v in sorted(counter.items())])
                cutoff = int(len(ds) * 0.7)
                if cutoff < len(ds) * 0.7:
                    cutoff += 1
                ps = [] # T-C mismatch
                for k, v in counter.items():
                    if v >= cutoff:
                        ps.append(k)
                
                is_nascent = False
                if len(ds) == 1:
                    is_nascent = len(ps) >= 5 and ds[0].tc_ratio >= 0.05
                elif len(ds) == 2:
                    is_nascent = len(ps) >= 3
                else:
                    is_nascent = len(ps) >= 2
                xt = "Nascent" if is_nascent else "Exists"
                for obj in ds:
                    obj.segment.set_tag("XT", xt)
                    
                t_count = None
                for a in ds:
                    if not a.segment.is_duplicate:
                        t_count = a.t_count
                assert t_count is not None
                
                # Metrics
                # BED6
                # 
                
                line = "\t".join(map(str, [
                    chrom, ds[0].start, ds[0].end, name, len(ds), ds[0].strand, 
                    xt, t_count, len(ps), s]))
                fw2.write(line + "\n")
                
            for obj in bundle.data:
                fw.write(obj.segment)
    f.close()
    fw.close()
    fw2.close()
    
    print("Finished!")
        
        
if __name__ == '__main__':
    main()
    
    