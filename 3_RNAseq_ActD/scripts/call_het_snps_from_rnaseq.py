#!/usr/bin/env python
import sys
from collections import defaultdict
import multiprocessing as mp
import numpy as np
from pyBioInfo.IO.File import BamFile
from pyBioInfo.Utils import BundleBuilder, SegmentTools


def call_het_snps(bamfile, chrom):
    # lines = []  
    snps = []  
    with BamFile(bamfile, random=True) as f:
        for bundle in BundleBuilder(f.fetch(chrom), keep=True):
            coverages = np.zeros(bundle.end_max - bundle.start_min, dtype=int)
            for align in bundle.data:
                for block in align.blocks:
                    for idx in range(block[0] - bundle.start_min, block[1] - bundle.start_min):
                        coverages[idx] += 1
            refs = dict()
            events = dict()
            for align in bundle.data:
                for e in SegmentTools.get_events(align.segment):
                    pos, ref, alt, qua, dis = e
                    if ref == "-": # insertion
                        continue
                    for idx1, ref1 in enumerate(ref):
                        pos1 = pos + idx1
                        if pos1 not in events:
                            refs[pos1] = ref1
                            events[pos1] = defaultdict(int)
                        events[pos1][alt] += 1
            for pos, counter in sorted(events.items()):
                cov = coverages[pos - bundle.start_min]
                ref = refs[pos]
                counter[ref] = cov - sum(counter.values())
                items = list(sorted(counter.items(), key=lambda item: item[1]))
                a1, c1 = items[-1]
                a2, c2 = items[-2]
                if a1 == "-" or a2 == "-":
                    continue
                r1, r2 = c1 / cov, c2 / cov
                b1 = (c1 >= 10 and r1 >= 0.1) or (c1 >= 5 and r1 >= 0.3)
                b2 = (c2 >= 10 and r2 >= 0.1) or (c2 >= 5 and r2 >= 0.3)
                if b1 and b2:
                    snp = [bundle.chrom, pos, ref, a1, a2]
                    snps.append(snp)
    return snps


def main():
    bamfile, sample, threads, outfile = sys.argv[1:]
    threads = int(threads)
    
    with BamFile(bamfile, random=True) as f:
        references = f.references

    results = []
    pool = mp.Pool(threads)
    for chrom in sorted(references):
        r = pool.apply_async(call_het_snps, (bamfile, chrom))
        results.append(r)
    pool.close()
    pool.join()
        
    with open(outfile, "w+") as fw:
        fw.write("##fileformat=VCFv4.2\n")
        for chrom in sorted(references):
            length = references[chrom]
            fw.write("##contig=<ID=%s,length=%d>\n" % (chrom, length))
        fw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        fw.write("##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase set for GT\">\n")
        fw.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample)
        for r in results:
            snps = r.get()
            for snp in snps:
                chrom, pos, ref, a1, a2 = snp
                alleles = [ref]
                if a1 not in alleles:
                    alleles.append(a1)
                if a2 not in alleles:
                    alleles.append(a2)
                alts = ",".join(alleles[1:])
                gt1 = alleles.index(a1)
                gt2 = alleles.index(a2)
                if gt1 > gt2:
                    gt1, gt2 = gt2, gt1
                gt = "%d|%d" % (gt1, gt2)
                line = "\t".join(map(str, [chrom, pos + 1, ".", ref, alts, ".", "PASS", ".", "GT", gt]))
                fw.write(line + "\n")
                    
    
if __name__ == "__main__":
    main()