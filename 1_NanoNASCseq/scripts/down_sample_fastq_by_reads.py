#!/usr/bin/env python
import optparse
import gzip
import numpy as np


def load_fastq(path):
    f = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    lines = list()
    for i, line in enumerate(f):
        lines.append(line)
        if i % 4 == 3:
            yield lines
            lines = list()
    f.close()
    
    
def down_sample_by_reads(infile, outfile, seed, num):
    np.random.seed(seed)
    reads1 = list(load_fastq(infile))
    idxs1 = np.arange(len(reads1))
    idxs2 = np.random.choice(idxs1, num, False)
    idxs2.sort()
    reads2 = [reads1[i] for i in idxs2]
    with open(outfile, "w+") as fw:
        for lines in reads2:
            for line in lines:
                fw.write(line)
    

def main():
    parser = optparse.OptionParser()
    parser.add_option("-i", "--infile", dest="infile", metavar="PATH", 
                      help="Input FASTQ file. [%default]")
    parser.add_option("-o", "--outfile", dest="outfile", metavar="PATH",
                      help="Output FASTQ file. [%default]")
    parser.add_option("-s", "--seed", dest="seed", metavar="INT", type="int", default=1, 
                      help="Random seed. [%default]")
    parser.add_option("-r", "--reads", dest="reads", metavar="INT", type="int", 
                      help="Number of reads at target FASTQ file. [%default]")
    options, _ = parser.parse_args()
    
    down_sample_by_reads(
        infile=options.infile,
        outfile=options.outfile,
        seed=options.seed,
        num=options.reads,
    )


if __name__ == "__main__":
    main()