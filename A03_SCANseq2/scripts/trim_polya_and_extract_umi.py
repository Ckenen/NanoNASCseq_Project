#!/usr/bin/env python
import sys
import numpy as np
import pysam

def find_polya(seq):
    m = np.zeros(len(seq) + 1)
    for i, base in enumerate(seq):
        if base == "A":
            m[i + 1] = m[i] + 1
        else:
            m[i + 1] = max(m[i] - 2, 0)
    vmax = max(m)
    start = 0
    end = None
    for i, v in enumerate(m):
        if v == 0:
            start = i
        if v == vmax:
            end = i
            break
    return start, end

def main():
    infile, outfile = sys.argv[1:]
    
    total = 0
    polya = 0
    with pysam.FastxFile(infile) as f, open(outfile, "w+") as fw:
        for i, read in enumerate(f):
            # if i >= 10000:
            #     break
            total += 1
            
            seq = read.sequence
            qua = read.quality
            
            w = min(100, len(seq))
            offset = len(seq) - w
            seq_tail = seq[-w:]
            x, y = find_polya(seq_tail)
            x, y = x + offset, y + offset
            
            a = y - x
            left = len(seq) - y
            if a >= 15 and left >= 6 and left <= 10:
                umi = seq[-8:]
                seq = seq[:x]
                qua = qua[:x]
                fw.write("@%s_%s\n%s\n+\n%s\n" % (read.name, umi, seq, qua))
            else:
                continue
            polya += 1
            
    print("Total reads:", total)
    print("PolyA reads:", polya)
    print("PolyA ratio:", polya / total)
    
    
if __name__ == "__main__":
    main()
        