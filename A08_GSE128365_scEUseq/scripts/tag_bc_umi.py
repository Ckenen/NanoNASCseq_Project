#!/usr/bin/env python
import sys
import gzip


def load_fastq(fqfile):
    with gzip.open(fqfile, "rt") as f:
        for i, line in enumerate(f):
            j = i % 4
            if j == 0:
                name = line.strip().split()[0][1:]
            elif j == 1:
                seq = line[:-1]
            elif j == 3:
                qua = line[:-1]
                yield name, seq, qua
            
def main():
    fqfile1, fqfile2 = sys.argv[1:]
    
    for read1, read2 in zip(load_fastq(fqfile1), load_fastq(fqfile2)):
        name1, seq1, qua1 = read1
        name2, seq2, qua2 = read2
        umi = seq1[0:6]
        bc = seq1[6:14]
        if len(seq1) < 40:
            continue        
        polyt = seq1[14:34]
        if polyt.count("T") < 19:
            continue
        print("@%s_%s_%s" % (bc, umi, name2))
        print(seq2)
        print("+")
        print(qua2)
    
    
if __name__ == "__main__":
    main()

