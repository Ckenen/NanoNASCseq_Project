#!/usr/bin/env python
import sys
from pyBioInfo.IO.File import BedFile


def main():
    with BedFile(sys.argv[1]) as f:
        for x in f:
            s = 'gene_id "%s"; transcript_id "%s";' % (x.name, x.name)
            row = [x.chrom, ".", "transcript", x.start + 1, x.end, ".", x.strand, ".", s]
            print("\t".join(map(str, row)))
            for block in x.blocks:
                row = [x.chrom, ".", "exon", block[0] + 1, block[1], ".", x.strand, ".", s]
                print("\t".join(map(str, row)))


if __name__ == "__main__":
    main()
