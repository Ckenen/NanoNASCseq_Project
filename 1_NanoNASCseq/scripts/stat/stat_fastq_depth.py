#!/usr/bin/env python
import sys
import subprocess


def load_fastq(path):
    if path.endswith(".gz"):
        p = subprocess.Popen(["gzip", "-dc", path], text=True, stdout=subprocess.PIPE)
        f = p.stdout
    else:
        f = open(path, "r")
    for i, line in enumerate(f):
        if i % 4 == 1:
            yield line.strip()
    f.close()
    

def main():
    infile, outfile = sys.argv[1:]
    
    with open(outfile, "w") as fw:
        fw.write("File\tReads\tBases\n")
        reads, bases = 0, 0
        for seq in load_fastq(infile):
            reads += 1
            bases += len(seq)
        fw.write(f"{infile}\t{reads}\t{bases}\n")

if __name__ == "__main__":    
    main()
    