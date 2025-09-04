#!/usr/bin/env python
import sys
from collections import Counter, defaultdict
import pandas as pd
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder, BedFile
from pyBioInfo.Utils import ShiftLoader, BlockTools


def is_fsm(ref, que):
    assert len(ref.blocks) > 1
    assert len(ref.blocks) > 1
    return ref.introns == que.introns


def is_ism(ref, que):
    assert len(ref.blocks) > 1
    assert len(ref.blocks) > 1
    if len(ref.introns) > len(que.introns):
        for i in range(len(ref.introns)):
            intron1 = ref.introns[i]
            intron2 = que.introns[0]
            if intron1 == intron2:
                if len(ref.introns) - i >= len(que.introns):
                    identical = True
                    for j in range(len(que.introns)):
                        if ref.introns[i + j] != que.introns[j]:
                            identical = False
                            break
                    if identical:
                        # check boundary
                        i1 = i
                        i2 = i + len(que.introns)
                        if i1 > 0:
                            if que.start - ref.introns[i1 - 1][1] < -10:
                                return False
                        if i2 < len(ref.introns):
                            if que.end - ref.introns[i2][0] >= 10:
                                return False
                        return True
                else:
                    return False
                break
    else:
        return False


def is_mono_exon_fsm(ref, que):
    assert len(ref.blocks) == 1
    assert len(ref.blocks) == 1
    start = max(ref.start, que.start)
    end = min(ref.end, que.end)
    length = end - start
    assert length > 0
    r = length / min(len(ref), len(que))
    return r >= 0.8


def is_mono_exon_ism(ref, que):
    assert len(ref.blocks) > 1
    assert len(que.blocks) == 1
    for block_start, block_end in ref.blocks:
        if min(block_end, que.end) - max(block_start, que.start) > 0:
            if que.start - block_start >= -10 and que.end - block_end <= 10:
                return True
            else:
                return False
    return False


def main():
    f_gtf, f_bed, f_tsv = sys.argv[1:]

    refs = []
    with GtfFile(f_gtf) as f:
        records = [x for x in f]
        for t in sorted(GtfTranscriptBuilder(records)):
            t.introns = tuple(BlockTools.gaps(t.blocks))
            t.gene_id = t.records["transcript"][0].attributes["gene_id"]
            t.transcript_id = t.records["transcript"][0].attributes["transcript_id"]
            refs.append(t)

    ques = []
    with BedFile(f_bed) as f:
        for umi in f:
            umi.introns = tuple(BlockTools.gaps(umi.blocks))
            ques.append(umi)

    loader = ShiftLoader(refs)
    rows = []
    for q in ques:
        for r in loader.fetch(obj=q):
            if r.strand == q.strand:
                if len(q.blocks) == 1:
                    if len(r.blocks) == 1:
                        if is_mono_exon_fsm(r, q):
                            rows.append([q.name, r.gene_id, r.transcript_id, "MonoExonFSM"])
                    else:
                        if is_mono_exon_ism(r, q):
                            rows.append([q.name, r.gene_id, r.transcript_id, "MonoExonISM"])
                else:
                    if len(r.blocks) == 1:
                        pass
                    else:
                        if is_fsm(r, q):
                            rows.append([q.name, r.gene_id, r.transcript_id, "FSM"])
                        elif is_ism(r, q):
                            rows.append([q.name, r.gene_id, r.transcript_id, "ISM"])
                        else:
                            pass

    dat = pd.DataFrame(rows, columns=["Name", "GeneID", "TranscriptID", "Category"])
    dat.to_csv(f_tsv, sep="\t", index=False)

    counter = defaultdict(int)
    for umi, d in dat.groupby(by="Name"):
        c = Counter(d["Category"])
        array = []
        for s in ["FSM", "ISM", "MonoExonFSM", "MonoExonISM"]:
            if c[s] == 1:
                array.append(s)
            elif c[s] > 1:
                array.append(s + "s")
        array.sort()
        counter["+".join(array)] += 1
    print("Combination\tCount")
    for k, v in sorted(counter.items()):
        print(k, v, sep="\t")
        
    
if __name__ == "__main__":
    main()