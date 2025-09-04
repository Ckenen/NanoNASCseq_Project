#!/usr/bin/env python
import sys
import os
import gzip
import glob
import subprocess


def main():
    infile, outdir = sys.argv[1:]

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    with open(infile) as f:
        fw = None
        for line in f:
            if line.startswith(">"):
                if fw is not None:
                    fw.close()
                chrom = line.strip("\n").split()[0][1:]
                fw = open(os.path.join(outdir, "%s.fasta" % chrom), "w+")
            fw.write(line)
        if fw is not None:
            fw.close()
            
    if True: # samtools faidx
        for path in glob.glob(os.path.join(outdir, "*.fasta")):
            cmd = "samtools faidx %s" % path
            subprocess.check_call(cmd, shell=True)


if __name__ == "__main__":
    main()