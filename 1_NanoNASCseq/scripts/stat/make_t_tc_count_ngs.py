#!/usr/bin/env python
import sys, pysam, optparse
from collections import defaultdict
import pandas as pd


MAPPER = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def load_segments(inbam, ignore_duplicate=True):
    with pysam.AlignmentFile(inbam) as f:
        for s in f:
            if s.is_duplicate and ignore_duplicate:
                continue
            yield s
        
        
def get_strand(segment):
    if segment.has_tag("ST"):
        return segment.get_tag("ST")
    else:
        return "-" if segment.is_reverse else "+"
     
                
def get_t(segment, strand):
    counter = defaultdict(int)
    for x in segment.get_tag("RC").split(";"):
        if x.strip() != "":
            k, v = x.split(",")
            if strand == "-":
                k = MAPPER[k]
            counter[k] = int(v)
    return counter["T"]


def get_tc(segment, strand):
    counter = defaultdict(int)
    for x in segment.get_tag("CE").split(";"):
        if x.strip() == "":
            continue
        pos, ref, alt, qua, dis = x.split(",")
        if strand == "-":
            ref, alt = MAPPER[ref], MAPPER[alt]
        counter["%s%s" % (ref, alt)] += 1
    return counter["TC"]


def make_t_tc_matrix(inbam, outfile):
    data = defaultdict(int)
    for segment in load_segments(inbam):
        strand = get_strand(segment)                
        t = get_t(segment, strand)
        tc = get_tc(segment, strand)
        data[(t, tc)] += 1
    rows = [[t, tc, count] for (t, tc), count in sorted(data.items())]
    m = pd.DataFrame(rows, columns=["#T", "TC", "Count"])
    m.to_csv(outfile, sep="\t" if outfile.endswith(".tsv") else ",", index=False)
            

def make_t_tc_matrix_main():
    parser = optparse.OptionParser(usage="%prog -i in.bam -o matrix.csv")
    parser.add_option("-i", "--ibam", dest="inbam")
    parser.add_option("-o", "--outfile", dest="outfile")
    options, _ = parser.parse_args()
    make_t_tc_matrix(
        inbam=options.inbam,
        outfile=options.outfile,
    )

    
if __name__ == "__main__":
    make_t_tc_matrix_main()