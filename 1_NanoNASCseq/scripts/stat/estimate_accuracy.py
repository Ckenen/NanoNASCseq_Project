#!/usr/bin/env python
import os, shutil, optparse
from collections import defaultdict
from Bio.Seq import Seq
import pysam
from spoa import poa
import pandas as pd
import subprocess
import multiprocessing as mp


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


def get_1_read_accuracy(segment):
    nm = segment.get_tag("NM")
    length = segment.infer_query_length()
    accuracy = 1 - nm / length
    return nm, length, accuracy


def get_2_read_accuracy(segment1, segment2):
    if segment1.get_tag("de") < segment2.get_tag("de"):
        segment = segment1
    else:
        segment = segment2
    return get_1_read_accuracy(segment)


def get_polished_sequence(name, seqs, tmpdir):
    fa1 = os.path.join(tmpdir, "%s.raw.fasta" % name)
    fa2 = os.path.join(tmpdir, "%s.consensus.fasta" % name)
    fa3 = os.path.join(tmpdir, "%s.polished.fasta" % name)
    sam = os.path.join(tmpdir, "%s.aligned.sam" % name)
    
    # make raw.fasta
    with open(fa1, "w+") as fw:
        for i, seq in enumerate(seqs):
            fw.write(">read%s\n" % i)
            fw.write("%s\n" % seq)
    
    # make consensus.fasta
    consensus, msa = poa(seqs)
    with open(fa2, "w+") as fw:
        fw.write(">%s\n" % name)
        fw.write("%s\n" % consensus)

    # make aligned.sam (minimap2 v2.28)
    cmd = "minimap2 -ax map-ont -o %s %s %s" % (sam, fa2, fa1)
    subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # make polished.fasta (racon v1.4.20)
    cmd = "racon %s %s %s > %s" % (fa1, sam, fa2, fa3)
    subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    polished_seq = None
    with pysam.FastxFile(fa3) as f:
        for record in f:
            polished_seq = record.sequence
            break

    os.remove(fa1)
    os.remove(fa2)
    os.remove(fa3)
    os.remove(sam)
    
    return name, polished_seq


def estimate_accuracy(bamfile, mmi, bedfile, threads, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    tmpdir = os.path.join(outdir, "tmp")
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
        
    rows = []
    tasks = []
    sizes = dict()
    with pysam.AlignmentFile(bamfile) as f:
        for chrom in f.references:
            segments = defaultdict(list)
            for s in f.fetch(chrom):
                segments[s.get_tag("CN")].append(s)
            for name, ss in segments.items():
                sizes[name] = len(ss)
                if len(ss) == 1:
                    nm, length, accuracy = get_1_read_accuracy(ss[0])
                    rows.append([name, len(ss), nm, length, accuracy])
                elif len(ss) == 2:
                    nm, length, accuracy = get_2_read_accuracy(ss[0], ss[1])
                    rows.append([name, len(ss), nm, length, accuracy])
                else:
                    seqs = []
                    for s in ss:
                        seq = s.query_sequence
                        if s.is_reverse:
                            seq = reverse_complement(seq)
                        seqs.append(seq)
                    tasks.append([name, seqs])

    results = []
    pool = mp.Pool(threads)
    for name, seqs in tasks:
        results.append(pool.apply_async(get_polished_sequence, (name, seqs, tmpdir)))
        # break
    pool.close()
    pool.join()
    
    polished_fa = os.path.join(outdir, "all.polished.fasta")
    with open(polished_fa, "w+") as fw:                
        for r in results:
            name, polished_seq = r.get()
            fw.write(">%s\n" % name)
            fw.write("%s\n" % polished_seq)
            
    polished_sam = os.path.join(outdir, "all.aligned.sam")
    cmd = "minimap2 -ax splice -u f -Y --MD -o %s --junc-bed %s -t %d %s %s" % (polished_sam, bedfile, threads, mmi, polished_fa)
    subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    with pysam.AlignmentFile(polished_sam) as f:
        for s in f:
            if s.is_unmapped or s.is_supplementary or s.is_secondary or s.mapping_quality < 60:
                continue
            name = s.query_name
            nm, length, accuracy = get_1_read_accuracy(s)
            rows.append([name, sizes[name], nm, length, accuracy])          
        
    m = pd.DataFrame(rows)
    m.columns = ["UMI", "Size", "NM", "Length", "Accuracy"]
    m.to_csv(os.path.join(outdir, "summary.csv"), index=False)
    
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
        
    
def main():
    parser = optparse.OptionParser(usage="%prog [options]")
    parser.add_option("-b", "--bam", dest="bamfile", help="Input BAM file. [%default]")
    parser.add_option("-m", "--mm2-index", dest="mmi", help="Minimap2 index. [%default]")
    parser.add_option("-g", "--gene", dest="bedfile", help="Gene annotation in BED format. [%default]")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, help="Threads. [%default]")
    parser.add_option("-o", "--outdir", dest="outdir", help="Output directory. [%default]")
    options, _ = parser.parse_args()
    estimate_accuracy(
        bamfile=options.bamfile,
        mmi=options.mmi,
        bedfile=options.bedfile,
        threads=options.threads,
        outdir=options.outdir,
    )


if __name__ == "__main__":
    main()