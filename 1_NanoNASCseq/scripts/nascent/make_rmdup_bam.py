#!/usr/bin/env python
import optparse
import pysam

# For NanoNASC-seq

def main():
    parser = optparse.OptionParser(usage="%prog [options] input.bam output.bam")
    parser.add_option("-r", "--min-reads", dest="min_reads", type="int", default=1, metavar="INT", 
                       help="Minimum reads of duplicate set. [%default]")
    options, args = parser.parse_args()
    f_bam, f_out = args
    min_reads = options.min_reads
    
    with pysam.AlignmentFile(f_bam) as inbam, \
        pysam.AlignmentFile(f_out, "wb", inbam) as outbam:
        for s in inbam:
            if s.get_tag("CS") < min_reads:
                continue
            if s.is_duplicate:
                continue
            outbam.write(s)
            
    
if __name__ == "__main__":
    main()
        