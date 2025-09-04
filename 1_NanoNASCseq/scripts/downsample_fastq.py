#!/usr/bin/env python
import sys
import gzip
import random
import optparse
import numpy as np


def load_reads(path):
    if path.endswith(".gz"):
        f = gzip.open(path, "rt")
    else:
        f = open(path)
    for i, line in enumerate(f):
        j = i % 4
        if j == 0:
            name = line.strip().split()[0][1:]
        elif j == 1:
            seq = line.strip()
        elif j == 3:
            qua = line.strip()
            yield name, seq, qua            
    f.close()
    

def main():
    parser = optparse.OptionParser(usage="%prog [options]")
    parser.add_option("-s", "--seed", dest="seed", type="int", default=1, metavar="INT", help="default: %default")
    parser.add_option("-r", "--source-reads", dest="source_reads", type="int", metavar="INT", help="default: %default")
    parser.add_option("-b", "--source-bases", dest="source_bases", type="int", metavar="INT", help="default: %default")
    parser.add_option("-t", "--target-bases", dest="target_bases", type="int", metavar="INT", help="default: %default")
    parser.add_option("-i", "--input", dest="infile", type="str", metavar="PATH", help="default: %default")
    parser.add_option("-o", "--output", dest="outfile", type="str", metavar="PATH", help="default: %default")
    options, args = parser.parse_args()

    seed = options.seed
    source_reads = options.source_reads
    source_bases = options.source_bases
    target_bases = options.target_bases
    infile = options.infile
    outfile = options.outfile
    
    assert source_reads > 0
    assert source_bases > 0
    assert target_bases < source_bases
    p = target_bases / source_bases
    
    np.random.seed(seed)
    if outfile is None:
        fw = sys.stdout
    elif outfile.endswith(".gz"):
        fw = gzip.open(outfile, "wt")
    else:
        fw = open(outfile, "w+")
    for name, seq, qua in load_reads(infile):
        if np.random.random() < p:
            fw.write("@%s\n%s\n+\n%s\n" % (name, seq, qua))
    fw.close()
    
    
if __name__ == "__main__":
    main()