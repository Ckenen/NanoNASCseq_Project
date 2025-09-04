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

    if infile.endswith(".gz"):
        f = gzip.open(infile, "rt")
    else:
        f = open(infile)
    fws = dict()
    for line in f:
        if line.startswith("#"):
            continue
        chrom = line.split("\t")[0]
        if chrom not in fws:
            fws[chrom] = open(os.path.join(outdir, "%s.gtf" % chrom), "w+")
        fws[chrom].write(line)
    for chrom, fw in fws.items():
        fw.close()
    f.close()

    if True: # bgzip and tabix
        for path in glob.glob(os.path.join(outdir, "*.gtf")):
            cmd1 = "bgzip -c %s > %s.gz" % (path, path)
            cmd2 = "tabix -p gff %s.gz" % path
            subprocess.check_call(cmd1, shell=True)
            subprocess.check_call(cmd2, shell=True)


if __name__ == "__main__":
    main()