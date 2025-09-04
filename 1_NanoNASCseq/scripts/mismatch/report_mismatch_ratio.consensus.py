#!/usr/bin/env python
import sys
from collections import Counter, defaultdict
import numpy as np
import pysam
from pyBioInfo.IO.File import FastaFile, Alignment


def parse_events(s):
    events = []
    for item in s.split(";"):
        if item == "":
            continue
        e = item.split(",")
        e[0], e[4] = int(e[0]), int(e[4])
        if "/" in e[3]:
            e[3] = [int(x) for x in e[3].split("/")]
        else:
            e[3] = int(e[3])
        events.append(e)
    return events


def main():
    f_bam, f_fasta, f_tsv, f_tsv2 = sys.argv[1:]

    mapper = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

    data = []

    with pysam.AlignmentFile(f_bam) as f, FastaFile(f_fasta, random=True) as fasta:
        for chrom in f.references:

            alignments = defaultdict(list)
            for s in f.fetch(chrom):
                obj = Alignment(s)            
                obj.ce = parse_events(s.get_tag("CE"))
                alignments[s.get_tag("CN")].append(obj)

            for cn, v in alignments.items():
                size = len(v)
                start = min([x.start for x in v])
                end = max([x.end for x in v])
                strand = v[0].strand

                covs = np.zeros(end - start, dtype=int)
                for a in v:
                    for x, y in a.blocks:
                        for idx in range(x - start, y - start):
                            covs[idx] += 1

                blocks = []
                idx0 = 0
                last_cov = 0
                for idx, cov in enumerate(covs):
                    if cov == 0:
                        if last_cov == 0:
                            pass
                        else:
                            blocks.append([idx0, idx])
                    else:
                        if last_cov == 0:
                            idx0 = idx
                        else:
                            pass
                    last_cov = cov
                if covs[-1] > 0:
                    blocks.append([idx0, len(covs)])
                blocks = [[x + start, y + start] for x, y in blocks]
                
                # blocks = [[x, y] for x, y in v[0].blocks]
                # if start >= blocks[0][1] or end <= blocks[-1][0]:
                #     msg = "\t".join(map(str, [chrom, start, end, cn, strand]))
                #     raise RuntimeError(msg)
                # blocks[0][0] = start
                # blocks[-1][1] = end

                seqs = [fasta.fast_fetch(chrom=chrom, start=x, end=y).upper() for x, y in blocks]
                base_counter = Counter("".join(seqs))

                event_counter = defaultdict(int)
                for a in v:
                    for e in a.ce:
                        assert start <= e[0] < end
                        event_counter[(e[0], e[1], e[2])] += 1

                confident_events = defaultdict(int)
                for (pos, ref, alt), count in event_counter.items():
                    cov = covs[pos - start]
                    if (count >= 0.75 * size) or (cov >= 2 and count >= 0.75 * cov):
                        confident_events["%s-%s" % (ref, alt)] += 1

                if strand == "-":
                    base_counter_new = defaultdict(int)
                    for base, count in base_counter.items():
                        base_counter_new[mapper[base]] = count
                    base_counter = base_counter_new
                    confident_events_new = defaultdict(int)
                    for m, count in confident_events.items():
                        confident_events_new["%s-%s" % (mapper[m[0]], mapper[m[2]])] = count
                    confident_events = confident_events_new

                data.append([cn, len(v), base_counter, confident_events])
             
    
    ms = []
    for ref in "ACTG":
        for alt in "ACTG":
            if ref != alt:
                ms.append("%s-%s" % (ref, alt))
                
    all_base_counter = defaultdict(int)
    all_event_counter = defaultdict(int)
    for cn, size, base_counter, event_counter in data:
        if size < 2:
            continue
        for base, count in base_counter.items():
            all_base_counter[base] += count
        for event, count in event_counter.items():
            all_event_counter[event] += count

    with open(f_tsv, "w+") as fw:
        fw.write("Type\tRefBase\tAltBase\tBaseCount\tEventCount\tRatio\n")
        for m in ms:
            t = "%s%s" % (m[0], m[2])
            base = all_base_counter[m[0]]
            count = all_event_counter[m]
            ratio = np.divide(count, base)
            line = "\t".join(map(str, [t, m[0], m[2], base, count, ratio]))
            fw.write(line + "\n")
        
    with open(f_tsv2, "w+") as fw:
        bases = ["A", "C", "G", "T"]
        mtypes = ms
        header = ["Name", "Size"]
        header.extend(bases)
        header.extend(mtypes)
        fw.write("\t".join(header) + "\n")
        for cn, size, base_counter, event_counter in data:
            row = [cn, size]
            row.extend([base_counter[b] for b in bases])
            row.extend([event_counter[t] for t in mtypes])
            fw.write("\t".join(map(str, row)) + "\n")
        
        
if __name__ == "__main__":
    main()