#!/usr/bin/env python
from collections import defaultdict
import optparse
import multiprocessing as mp
from pyBioInfo.IO.File import BamFile, FastaFile, BedFile
from pyBioInfo.Utils import BlockTools


def get_motif(f, chrom, start, end, strand):
    seq1 = list(f.fetch(chrom, start, start + 2))[0].sequence
    seq2 = list(f.fetch(chrom, end - 2, end))[0].sequence
    if strand == "-":
        s1 = FastaFile.reverse_complement(seq1)
        s2 = FastaFile.reverse_complement(seq2)
        seq1, seq2 = s2, s1
    motif = "%s-%s" % (seq1, seq2)
    return motif


def load_known_introns(bedfile, chrom):
    introns = set()
    with BedFile(bedfile, random=True) as f:
        try:
            for t in f.fetch(chrom):
                for x, y in BlockTools.gaps(t.blocks):
                    introns.add((t.chrom, x, y, t.strand))
        except ValueError:
            pass
    return introns
                

def worker(bamfile, bedfile, fafile, chrom, direction):
    data = defaultdict(list)
    with BamFile(bamfile, random=True) as f:
        for a in f.fetch(chrom):
            s = a.segment
            introns = list(BlockTools.gaps(a.blocks))
            strand = a.strand
            if s.is_paired and s.is_read2:
                strand = "-" if strand == "+" else "+"
            if direction == "R":
                strand = "-" if strand == "+" else "+"
            for x, y in introns:
                data[(chrom, x, y, strand)].append(a.name)
        
    counter = dict()
    for k, names in data.items():
        counter[k] = len(set(names))

    motifs = dict()
    if len(counter) > 0 and fafile is not None:
        with FastaFile(fafile, random=True) as f:
            for k in data.keys():
                motifs[k] = get_motif(f, *k)
        
    known = dict()
    if len(counter) > 0 and bedfile is not None:
        introns = load_known_introns(bedfile, chrom)
        for k in data.keys():
            known[k] = k in introns
    
    return counter, motifs, known
        
            
def main():
    parser = optparse.OptionParser(usage="%prog [options] <input.bam> <out.bed>")
    parser.add_option("-b", "--bed", dest="bed", help="Transcript annotation.")
    parser.add_option("-f", "--fasta", dest="fasta", help="Genome sequences.")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=0)
    parser.add_option("-d", "--direction", dest="direction", default="F")
    options, args = parser.parse_args()
    bamfile, outfile = args
    bedfile = options.bed
    fafile = options.fasta
    threads = options.threads
    direction = options.direction

    results = []
    pool = mp.Pool(threads)
    with BamFile(bamfile, random=True) as f:
        chroms = list(sorted(f.references.keys()))
    for chrom in chroms:
        r = pool.apply_async(worker, (bamfile, bedfile, fafile, chrom, direction))
        results.append(r)
    pool.close()
    pool.join()

    with open(outfile, "w+") as fw:
        fw.write("Chrom\tStart\tEnd\tStrand\tMotif\tKnown\tCount\n")
        for r in results:
            counter, motifs, knowns = r.get()
            for k in sorted(counter.keys()):
                row = [*k, motifs.get(k), knowns.get(k), counter[k]]
                fw.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()