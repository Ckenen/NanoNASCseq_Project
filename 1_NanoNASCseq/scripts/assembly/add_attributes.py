#!/usr/bin/env python
import sys
from collections import defaultdict
import pandas as pd
from pyBioInfo.IO.File import GtfFile


def load_gtf(path):
    records = [x for x in GtfFile(path)]
    tid2records = defaultdict(list)
    for r in records:
        if r.feature == "gene":
            continue
        tid2records[r.attributes["transcript_id"]].append(r)
    return records, tid2records


def make_gene_info_table(records):
    rows = []
    for r in records:
        if r.feature != "gene":
            continue
        rows.append([
            r.attributes["gene_id"], 
            r.attributes["gene_name"], 
            r.attributes["gene_type"], 
            r.chrom, 
            r.start, 
            r.end, 
            r.strand])
    gene_info = pd.DataFrame(rows, columns=[
        "GeneID", "GeneName", "GeneType", 
        "Chrom", "Start", "End", "Strand"])
    return gene_info
    
    
def make_transcript_info_table(records, gene_info):
    gid2gname = dict()
    gid2gtype = dict()
    for gid, gname, gtype in gene_info[["GeneID", "GeneName", "GeneType"]].values:
        gid2gname[gid] = gname
        gid2gtype[gid] = gtype
        
    novel_isoform_counter = defaultdict(int)
    
    rows = []
    for r in records:
        if r.feature != "transcript":
            continue
        gid = r.attributes["gene_id"]
        gname = gid2gname[gid]
        gtype = gid2gtype[gid]
        tid = r.attributes["transcript_id"]
        tname = r.attributes.get("transcript_name")
        if tname is None:
            # The name of the first novel isoform is "gname-501"
            tname = "%s-%d" % (gname, 501 + novel_isoform_counter[gname])
            novel_isoform_counter[gname] += 1
        ttype = r.attributes.get("transcript_type")
        if ttype is None:
            ttype = "Unknown"
        rows.append([
            tid, tname, ttype,
            gid, gname, gtype,
            r.chrom, r.start, r.end, r.strand])
    transcript_info = pd.DataFrame(rows, columns=[
        "TranscriptID", "TranscriptName", "TranscriptType", 
        "GeneID", "GeneName", "GeneType", 
        "Chrom", "Start", "End", "Strand"])
    return transcript_info


def update_transcript_records(transcript_info, tid2records):
    columns = ["TranscriptID", "TranscriptName", "TranscriptType", "GeneName", "GeneType"]
    for tid, tname, ttype, gname, gtype in transcript_info[columns].values:
        for r in tid2records[tid]:
            r.attributes["gene_name"] = gname
            r.attributes["gene_type"] = gtype
            r.attributes["transcript_name"] = tname
            r.attributes["transcript_type"] = ttype
    

def main():
    infile, outfile = sys.argv[1:]
    
    records, tid2records = load_gtf(infile)
    
    gene_info = make_gene_info_table(records)
    gene_info.to_csv(outfile + ".gene_info.csv", index=False)
        
    transcript_info = make_transcript_info_table(records, gene_info)
    transcript_info.to_csv(outfile + ".transcript_info.csv", index=False)
        
    update_transcript_records(transcript_info, tid2records)
    
    with open(outfile, "w+") as fw:
        for r in records:
            fw.write(r.format() + "\n")
            
        
if __name__ == "__main__":
    main()