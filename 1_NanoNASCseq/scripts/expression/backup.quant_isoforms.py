#!/usr/bin/env python
import sys
import pickle
from collections import defaultdict
from pyBioInfo.IO.File import BedFile
from pyBioInfo.Utils import ShiftLoader, BlockTools


def compare(umi, transcript):
    if umi.strand != transcript.strand:
        return None
    if len(umi.introns) == 0:
        return None
    if len(umi.introns) == len(transcript.introns):
        if umi.introns == transcript.introns:
            return "FSM"
    elif len(umi.introns) < len(transcript.introns):
        n = len(umi.introns)
        if transcript.strand == "+":
            conflict = False
            for i in range(n):
                if umi.introns[-i-1] != transcript.introns[-i-1]:
                    conflict = True
                    break
            if not conflict:
                if umi.start - transcript.introns[-n-1][1] >= -10:
                    return "ISM"

        else:
            conflict = False
            for i in range(n):
                if umi.introns[i] != transcript.introns[i]:
                    conflict = True
                    break
            if not conflict:
                if transcript.introns[n][0] - umi.end >= -10:
                    return "ISM"
    return None
    

def main():
    f_anno, f_umi, f_info, f_parental, prefix = sys.argv[1:]

    umi_properties = dict()

    for line in open(f_info):
        row = line.strip("\n").split("\t")
        name, size, t_count, tc_count = row[3], int(row[4]), int(row[7]), int(row[8])
        tmp = {"Size": size, "T": t_count, "TC": tc_count}
        umi_properties[name] = tmp
        
    with open(f_parental) as f:
        for line in f:
            name, parental = line.strip("\n").split("\t")
            umi_properties[name]["Parental"] = parental 
            
    umis = []
    with BedFile(f_umi) as f:
        for umi in f:
            d = umi_properties[umi.name[4:]]
            umi.size = d["Size"]
            umi.t_count = d["T"]
            umi.tc_count = d["TC"]
            umi.parental = d["Parental"]
            umi.introns = tuple(BlockTools.gaps(umi.blocks))
            umis.append(umi)
                    
    with open(f_anno, "rb") as f:
        transcripts = pickle.load(f)

    transcripts = list(filter(lambda item: len(item.blocks) >= 2, transcripts))
    umis = list(filter(lambda item: item.size >= 2, umis))

    loader = ShiftLoader(transcripts)
    for umi in umis:
        ts = defaultdict(list)
        for t in loader.fetch(obj=umi):
            s = compare(umi, t)
            if s == "FSM" or s == "ISM":
                ts[s].append(t)
        umi.transcripts = ts
        
    counter = defaultdict(int)
    for umi in umis:
        ts1 = umi.transcripts["FSM"]
        ts2 = umi.transcripts["ISM"]
        counter[(len(ts1), len(ts2))] += 1
    
    with open(prefix + ".pkl", "wb") as fw:
        pickle.dump(counter, fw)
        
    counts = dict()
    for t in transcripts:
        counts[t.name] = [0, 0, 0, 0, 0, 0]
    for umi in umis:
        if umi.size < 2:
            continue
        if len(umi.introns) == 0:
            continue
        ts1, ts2 = umi.transcripts["FSM"], umi.transcripts["ISM"]
        t = None
        if len(ts1) == 1 and len(ts2) == 0:
            t = ts1[0]
        if len(ts1) == 0 and len(ts2) == 1:
            t = ts2[0]
        if t is None:
            continue
        counts[t.name][0] += 1
        if umi.parental == "P":
            counts[t.name][1] += 1
        elif umi.parental == "M":
            counts[t.name][2] += 1
        if umi.tc_count >= 2:
            counts[t.name][3] += 1
            if umi.parental == "P":
                counts[t.name][4] += 1
            elif umi.parental == "M":
                counts[t.name][5] += 1
    
    with open(prefix + ".tsv", "w+") as fw:
        fw.write("TranscriptID\tGeneID\tGeneName\tTotal\tTotal.P\tTotal.M\tNascent\tNascent.P\tNascent.M\n")
        for t in transcripts:
            vs = counts[t.name]
            row = [t.name, t.gene_id, t.gene_name, vs[0], vs[1], vs[2], vs[3], vs[4], vs[5]]
            fw.write("\t".join(map(str, row)) + "\n")
            

if __name__ == "__main__":
    main()