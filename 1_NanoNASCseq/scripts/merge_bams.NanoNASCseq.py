#!/usr/bin/env python
import optparse, subprocess, pysam

USAGE="""
%prog -s 2 --rmdup filelist.txt merged.bam

"""


def run_cmd(cmd):
    subprocess.check_call(cmd, shell=True)    
    

def main():
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option("-r", "--rmdup", dest="rmdup", action="store_true", default=False)
    parser.add_option("-s", "--size", dest="read_size", type="int", default=1)
    options, args = parser.parse_args()
    rmdup = options.rmdup
    min_reads = options.read_size
    filelist, outfile = args
    
    paths = [line.strip("\n") for line in open(filelist)]
    
    tmpfile = outfile.replace(".bam", ".unsort.bam")
    statfile = outfile.replace(".bam", ".flagstat")
    
    fw = None
    for path in paths:
        with pysam.AlignmentFile(path) as f:
            if fw is None:
                fw = pysam.AlignmentFile(tmpfile, "wb", f)
            for s in f:
                if s.get_tag("CS") < min_reads:
                    continue
                if rmdup and s.is_duplicate:
                    continue
                fw.write(s)
    fw.close()
    
    run_cmd(f"samtools sort -@ 8 -o {outfile} {tmpfile}")
    run_cmd(f"samtools index -@ 8 {outfile}")
    run_cmd(f"rm {tmpfile}")
    run_cmd(f"samtools flagstat -@ 8 {outfile} > {statfile}")
    
    
if __name__ == "__main__":
    main()