#!/usr/bin/env python
import sys
from collections import defaultdict
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import Alignment
from pyBioInfo.Utils import ShiftLoader

def load_phased_snps(f, chrom):
    if isinstance(f, pysam.VariantFile):
        sample = list(f.header.samples)[0]
        try:
            for record in f.fetch(chrom):
                start = record.pos - 1
                ps = record.samples[sample]["PS"]
                if ps == ".":
                    continue
                if ps != "PATMAT" and ps != "0": # GIAB VCF
                    continue
                gt = record.samples[sample]["GT"]
                allele1 = record.alleles[gt[0]]
                allele2 = record.alleles[gt[1]]
                if len(allele1) > 1 or len(allele2) > 1:
                    continue
                if allele1 == allele2:
                    continue
                snp = GRange(chrom=chrom, start=start, end=start + 1)
                snp.allele1 = allele1
                snp.allele2 = allele2
                yield snp
        except ValueError:
            pass
    elif isinstance(f, pysam.TabixFile):
        try:
            for line in f.fetch(chrom):
                chrom, start, end, sample = line.strip("\n").split("\t")[:4]
                start, end = int(start), int(end)
                base_p, base_m = sample.split("|")
                if base_p == base_m or len(base_p) > 1 or len(base_m) > 1:
                    continue
                snp = GRange(chrom=chrom, start=start, end=end)
                snp.allele1 = base_p
                snp.allele2 = base_m
                yield snp
        except ValueError:
            pass
    else:
        raise RuntimeError()


def main():
    f_bam, f_vcf = sys.argv[1:]

    print("Name\tSize\tAlleles")

    with pysam.AlignmentFile(f_bam) as f, pysam.VariantFile(f_vcf) as f2:

        for chrom in f.references:
            data = defaultdict(list)
            for s in f.fetch(chrom):
                data[s.get_tag("CN")].append(Alignment(s))
            
            array = []
            for k, v in data.items():
                obj = GRange(
                    chrom=v[0].chrom, 
                    start=min([a.start for a in v]), 
                    end=max([a.end for a in v]), 
                    strand=v[0].strand, name=k)
                obj.alignments = v
                array.append(obj)
            array.sort()

            loader = ShiftLoader(load_phased_snps(f2, chrom))
            for x in array:
                snps = list(loader.fetch(chrom=x.chrom, start=x.start, end=x.end))
                alleles = []
                for snp in snps:
                    allele = None
                    counter = defaultdict(int)
                    for a in x.alignments:
                        if a.start <= snp.start < a.end:
                            base = a.get_query_base(snp.start)
                            if base is None:
                                continue
                            counter[base] += 1
                    c1 = counter[snp.allele1]
                    c2 = counter[snp.allele2]
                    t = sum(counter.values())
                    if t > 0:
                        if c1 > t * 0.5:
                            allele = "_".join(map(str, [snp.chrom, snp.start, snp.allele1, snp.allele2, snp.allele1]))
                        elif c2 > t * 0.5:
                            allele = "_".join(map(str, [snp.chrom, snp.start, snp.allele1, snp.allele2, snp.allele2]))
                    if allele is not None:
                        alleles.append(allele)
                print(x.name, len(x.alignments), ";".join(alleles), sep="\t")


if __name__ == "__main__":
    main()