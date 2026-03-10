#!/usr/bin/env python
import sys
import numpy as np
import pyBigWig


def main():
    bwfile1, bwfile2, width, outfile = sys.argv[1:]
    width = int(width)    
    step = int(width / 2)
    fw = open(outfile, 'w+')
    fw.write('BinID\tCoverage\n')
    for strand, path in zip(["+", "-"], [bwfile1, bwfile2]):
        f = pyBigWig.open(path)
        for chrom, length in f.chroms().items():
            n = int(length / step)
            if length % step > 0:
                n += 1
            covs = np.zeros(n)
            for item in f.intervals(chrom):
                pos = item[1]
                i = int(pos / step)
                covs[i] += item[2]
                if i > 0:
                    covs[i-1] += item[2]
            for i, cov in enumerate(covs):
                if cov == 0:
                    continue
                start = i * step
                end = min(i * step + width, length)
                fw.write('%s:%s-%s(%s)\t%s\n' % (chrom, start, end, strand, cov))
        f.close()
    fw.close()
            
    
if __name__ == "__main__":
    main()