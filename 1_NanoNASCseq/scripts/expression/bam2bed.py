#!/usr/bin/env python
import sys
import pysam
from pyBioInfo.IO.File import Alignment


def main():
    f_bam = sys.argv[1]    
    with pysam.AlignmentFile(f_bam) as f:
        for s in f:
            if s.is_duplicate:
                continue
            a = Alignment(s)
            a.name = s.get_tag("CN")
            print(a.format("BED"))


if __name__ == "__main__":
    main()