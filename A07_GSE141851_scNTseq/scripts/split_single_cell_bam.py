#!/usr/bin/env python
import sys
import os
import subprocess
from collections import defaultdict
import pysam


def main():
    inbam, outdir = sys.argv[1:]
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    with pysam.AlignmentFile(inbam) as f:
        counter = defaultdict(int)
        for s in f:
            counter[s.get_tag("XC")] += 1

    if True:
        num_cells = 800
        items = list(sorted(counter.items(), key=lambda x: x[1], reverse=True))
        barcodes = [item[0] for item in items[:num_cells]]
    else:
        min_reads = 2000
        barcodes = set([item[0] for item in filter(lambda item: item[1] >= min_reads, counter.items())])
    print("Number of cells:", len(barcodes))

    with open(outdir + "/barcodes.txt", "w+") as fw:
        for bc in sorted(barcodes):
            fw.write("%s\t%d\n" % (bc, counter[bc]))

    with pysam.AlignmentFile(inbam) as f:
        handles = dict()
        for s in f:
            barcode = s.get_tag("XC")
            if barcode in barcodes:
                if barcode not in handles:
                    handles[barcode] = pysam.AlignmentFile(outdir + "/%s.unsorted.bam" % barcode, "wb", f)
                handles[barcode].write(s)
        for h in handles.values():
            h.close()

    for barcode in barcodes:
        path1 = outdir + "/%s.unsorted.bam" % barcode
        path2 = outdir + "/%s.bam" % barcode
        cmd = "samtools sort -o %s %s" % (path2, path1)
        subprocess.check_call(cmd, shell=True)
        cmd = "samtools index %s" % path2
        subprocess.check_call(cmd, shell=True)
        cmd = "rm %s" % path1
        subprocess.check_call(cmd, shell=True)
        
    
if __name__ == "__main__":
    main()