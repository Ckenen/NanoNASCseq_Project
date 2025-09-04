#!/usr/bin/env python
import sys
import pysam


def main():
    infile, outfile = sys.argv[1:]
    
    with pysam.AlignmentFile(infile) as f, pysam.AlignmentFile(outfile, "wb", f) as fw:
        for s in f.fetch(until_eof=True):
            name = s.query_name
            umi = name.split(":")[-1]
            name = name[:-(len(umi) + 1)]
            s.query_name = name
            s.set_tag("UM", umi)
            fw.write(s)
    
    
if __name__ == "__main__":
    main()