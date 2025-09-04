#!/usr/bin/env python
import sys
import os
from collections import defaultdict
import edlib
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from pyBioInfo.IO.File import FastqFile
from pygz import PigzFile


# debug
DEBUG = False

     
def get_cell(fqfile):
    cell = fqfile.split("/")[-1]
    if cell.endswith(".gz"):
        cell = cell[:-3]
    if cell.endswith(".fastq"):
        cell = cell[:-6]
    if cell.endswith(".fq"):
        cell = cell[:-3]
    return cell


def load_fastq(path):
    with FastqFile(path) as f:
        for r in f:
            yield Read(r.name, r.sequence, r.quality)


def edlib_align(que, ref):
    a = edlib.align(que, ref, task="locations", mode="HW")
    x, y = a["locations"][0]
    return x, y + 1, a["editDistance"]


def write_fastq(fw, name, seq, qua):
    fw.write("@%s\n%s\n+\n%s\n" % (name, seq, qua))


def calculate_percentage(n1, n2):
    return 0.0 if n2 == 0 else n1 * 100 / n2


def get_ratio(n1, n2):
    return 0.0 if n2 == 0 else n1 / n2


def reverse_complement(seq, qua):
    return str(Seq(seq).reverse_complement()), qua[::-1]
    
      
def find_polya(seq): # seq is partial sequence
    if len(seq) == 0:
        return 0, 0
    scores = np.zeros(len(seq) + 1, dtype=int)
    for i, base in enumerate(seq):
        if base == "A":
            scores[i + 1] = scores[i] + 1
        else:
            scores[i + 1] = max(scores[i] - 2, 0)
    vmax = np.max(scores)
    x, y = 0, 0
    for i, score in enumerate(scores):
        if score == 0:
            x = i
        if score == vmax:
            y = i
            break
    return x, y


def find_full_polya(seq): # seq = mRNA + polyA
    assert len(seq) > 0
    score = 0
    vmax = -10
    idx = -1
    for i, base in enumerate(seq[::-1]):
        if base == "A":
            score += 1
        else:
            score -= 2
        if score > vmax:
            vmax = score
            idx = i
        elif score < -10:
            break
        elif score < vmax - 10:
            break
    x = len(seq) - idx - 1
    y = len(seq)
    return x, y


class Read(object):
    def __init__(self, name, seq, qua):
        self.name = name
        self.raw_seq = seq
        self.raw_qua = qua
        self.seq = seq # trimmed seq
        self.qua = qua
        self.raw_length = len(self.raw_seq)
        self.trim_length = None

        self.min_raw_length = 200
        self.min_trim_length = 200
                
        self.status = True
        self.reason = None
        
        self.tso_seq_f = "AAGCAGTGGTATCAACGCAGAGTAC"
        self.tso_seq_r = "GTACTCTGCGTTGATACCACTGCTT"
        self.expected_anchor_seq = "TTGCATCT"
        self.expected_anchor_umi_length = 20
        
        self.tso_h_ed = None
        self.tso_t_ed = None
        self.n_polyt = None
        self.n_polya = None
        self.direction = None
        self.anchor_umi_seq = None
        self.anchor_x = None
        self.anchor_ed = None
        self.anchor_seq = None
        self.umi_seq = None
        
    def _check_raw_length(self):
        if self.raw_length < self.min_raw_length:
            self.status = False
            self.reason = "RawTooShort"
    
    def _find_and_trim_tso(self):
        w = 30 # int(len(tso_seq_f) * 1.2)
        offset = len(self.seq) - w
        x1, y1, ed1 = edlib_align(self.tso_seq_f, self.seq[:w])
        x2, y2, ed2 = edlib_align(self.tso_seq_r, self.seq[-w:])
        x2, y2 = x2 + offset, y2 + offset
        assert x2 > y1
        if ed1 <= 8 and ed2 <= 8: # ED of TSO
            self.seq = self.seq[y1:x2]
            self.qua = self.qua[y1:x2]
        else:
            self.status = False
            self.reason = "NoTSO"
        self.tso_h_ed = ed1
        self.tso_t_ed = ed2        
    
    def _check_chimeric(self):
        # x1, y1, ed1 = edlib_align("AAGCAGTGGTATCAACGCAGAGTACATGGG", seq)
        # x2, y2, ed2 = edlib_align("CCCATGTACTCTGCGTTGATACCACTGCTT", seq)
        x1, y1, ed1 = edlib_align(self.tso_seq_f, self.seq)
        x2, y2, ed2 = edlib_align(self.tso_seq_r, self.seq)
        if ed1 <= 6 or ed2 <= 6: # ED of chimeric TSO
            self.status = False
            self.reason = "IsChimeric"
            
    def _determine_direction(self):
        w = self.expected_anchor_umi_length + 40
        offset = len(self.seq) - w
        seq_f, qua_f = self.seq, self.qua
        seq_r, qua_r = reverse_complement(self.seq, self.qua)
        x1, y1 = find_polya(seq_f[-w:])
        x2, y2 = find_polya(seq_r[-w:])
        x1, y1 = x1 + offset, y1 + offset
        x2, y2 = x2 + offset, y2 + offset
        self.n_polya, self.n_polyt = y1 - x1, y2 - x2
        
        self.direction = "U"
        if self.n_polya >= 15 and self.n_polyt < 10: # PolyA length
            seq, qua = seq_f, qua_f
            x3, y3 = x1, y1
            self.direction = "F"
        elif self.n_polyt >= 15 and self.n_polya < 10:
            seq, qua = seq_r, qua_r
            x3, y3 = x2, y2
            self.direction = "R"
            
        if self.direction != "U":
            seq3, qua3 = seq[:y3], qua[:y3] # mRNA + polyA, remove ATGGG and polyA
            x4, y4 = find_full_polya(seq3)
            self.seq = seq3[5:x4] # mRNA
            self.qua = qua3[5:x4]
            self.trim_length = len(self.seq)
            # self.seq = seq3
            # self.qua = qua3
            self.anchor_umi_seq = seq[y3:] # anchor + umi
        else:
            self.status = False
            self.reason = "NoDirection"
                
    def _check_anchor(self):
        if len(self.anchor_umi_seq) + 2 < len(self.expected_anchor_seq):
            self.status = False # too short
            self.reason = "NoAnchor"
        else:
            x, y, ed = edlib_align(self.expected_anchor_seq, self.anchor_umi_seq)
            if x <= 3 and ed <= 2: # Start of anchor <= 3, ED of anchor <= 2
                self.umi_seq = self.anchor_umi_seq[y:]
                self.anchor_seq = self.anchor_umi_seq[x:y]
            else:
                self.status = False
                self.reason = "NoAnchor"
            self.anchor_x = x
            self.anchor_ed = ed
    
    def _check_umi(self):
        if len(self.umi_seq) < 11 or len(self.umi_seq) > 13: # UMI length
            self.status = False
            self.reason = "NoUMI"
    
    def _check_trim_length(self):
        if self.trim_length < 200: # Final length
            self.status = False
            self.reason = "TrimTooShort"
    
    def trim(self):
        if self.status:
            self._check_raw_length()
        if self.status:
            self._find_and_trim_tso()
        if self.status:
            self._check_chimeric()
        if self.status:
            self._determine_direction()
        if self.status:
            self._check_anchor()
        if self.status:
            self._check_umi()
        if self.status:
            self._check_trim_length()
       
       
def main():
    infile, outdir = sys.argv[1:]
    
    cell = get_cell(infile)
    print("Infile: %s" % infile)
    print("Outfile: %s" % outdir)
    print("Cell: %s" % cell)
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    # Check cell information
    # if False:
    #     meta = pd.read_excel("data/NanoNASCseq_All.xlsx")
    #     assert meta[meta["Cell"] == cell]["LibStruct"].values[0] == "struct2"
    #     assert meta[meta["Cell"] == cell]["UMI"].values[0] == 20
            
    raw_length_counter = defaultdict(int)
    tso_ed_counter = defaultdict(int)
    polya_counter = defaultdict(int)
    anchor_x_ed_counter = defaultdict(int)
    umi_len_counter = defaultdict(int)
    trim_length_counter = defaultdict(int)
    reason_counter = defaultdict(int)
    n_fail = 0
    n_pass = 0
    with FastqFile(infile) as f, \
        open(os.path.join(outdir, "trimmed.fastq"), "w+") as fw:
        for i, r in enumerate(f):
            # if i >= 10000:
            #     break
            read = Read(r.name, r.sequence, r.quality)
            read.trim()
            if read.status:
                n_pass += 1
                write_fastq(fw, "%s:%s" % (read.name, read.umi_seq), read.seq, read.qua)
            else:
                reason_counter[read.reason] += 1
                n_fail += 1
            if read.raw_length is not None:
                raw_length_counter[read.raw_length] += 1
            if read.tso_h_ed is not None:
                tso_ed_counter[(read.tso_h_ed, read.tso_t_ed)] += 1
            if read.n_polya is not None:
                polya_counter[(read.n_polya, read.n_polyt)] += 1
            if read.anchor_x is not None:
                anchor_x_ed_counter[(read.anchor_x, read.anchor_ed)] += 1
            if read.umi_seq is not None:
                umi_len_counter[len(read.umi_seq)] += 1
            if read.trim_length is not None:
                trim_length_counter[read.trim_length] += 1
    
    # Report    
    n_total = n_pass + n_fail
    print("-" * 80)
    print("Total reads: %d (%.2f%%)" % (n_total, calculate_percentage(n_total, n_total)))
    print("Raw too short: %d (%.2f%%)" % (reason_counter["RawTooShort"], calculate_percentage(reason_counter["RawTooShort"], n_total)))
    print("No tso: %d (%.2f%%)" % (reason_counter["NoTSO"], calculate_percentage(reason_counter["NoTSO"], n_total)))
    print("Is chimeric: %d (%.2f%%)" % (reason_counter["IsChimeric"], calculate_percentage(reason_counter["IsChimeric"], n_total)))
    print("No direction: %d (%.2f%%)" % (reason_counter["NoDirection"], calculate_percentage(reason_counter["NoDirection"], n_total)))
    print("No anchor: %d (%.2f%%)" % (reason_counter["NoAnchor"], calculate_percentage(reason_counter["NoAnchor"], n_total)))
    print("No umi: %d (%.2f%%)" % (reason_counter["NoUMI"], calculate_percentage(reason_counter["NoUMI"], n_total)))
    print("Trim too short: %d (%.2f%%)" % (reason_counter["TrimTooShort"], calculate_percentage(reason_counter["TrimTooShort"], n_total)))
    print("Pass: %d (%.2f%%)" % (n_pass, calculate_percentage(n_pass, n_total)))
    
    with open(os.path.join(outdir, "raw_length.tsv"), "w+") as fw:
        fw.write("Length\tNumber\tRatio\n")
        t = sum(raw_length_counter.values())
        for k, n in sorted(raw_length_counter.items()):
            fw.write("\t".join(map(str, [k, n, get_ratio(n, t)])) + "\n")
            
    with open(os.path.join(outdir, "tso.tsv"), "w+") as fw:
        fw.write("ED1\tED2\tNumber\tRatio\n")
        t = sum(tso_ed_counter.values())
        for (ed1, ed2), n in sorted(tso_ed_counter.items()):
            fw.write("\t".join(map(str, [ed1, ed2, n, get_ratio(n, t)])) + "\n") 
            
    with open(os.path.join(outdir, "polya.tsv"), "w+") as fw:
        fw.write("PolyA\tPolyT\tNumber\tRatio\n")
        t = sum(polya_counter.values())
        for (n1, n2), n in sorted(polya_counter.items()):
            fw.write("\t".join(map(str, [n1, n2, n, get_ratio(n, t)])) + "\n")
            
    with open(os.path.join(outdir, "anchor.tsv"), "w+") as fw:
        fw.write("Start\tED\tNumber\tRatio\n")
        t = sum(anchor_x_ed_counter.values())
        for (x, ed), n in sorted(anchor_x_ed_counter.items()):
            fw.write("\t".join(map(str, [x, ed, n, get_ratio(n, t)])) + "\n")
    
    with open(os.path.join(outdir, "umi.tsv"), "w+") as fw:
        fw.write("Length\tNumber\tRatio\n")
        t = sum(umi_len_counter.values())
        for k, n in sorted(umi_len_counter.items()):
            fw.write("\t".join(map(str, [k, n, get_ratio(n, t)])) + "\n")
    
    with open(os.path.join(outdir, "trim_length.tsv"), "w+") as fw:
        fw.write("Length\tNumber\tRatio\n")
        t = sum(trim_length_counter.values())
        for k, n in sorted(trim_length_counter.items()):
            fw.write("\t".join(map(str, [k, n, get_ratio(n, t)])) + "\n")

    with open(os.path.join(outdir, "stats.tsv"), "w+") as fw:
        fw.write("File\tTotal\tRawTooShort\tNoTSO\tIsChimeric\tNoDirection\tNoAnchor\tNoUMI\tTrimTooShort\tPass\tPassRatio\n")
        fw.write("\t".join(map(str, [os.path.basename(infile), 
                                     n_total, 
                                     reason_counter["RawTooShort"], 
                                     reason_counter["NoTSO"], 
                                     reason_counter["IsChimeric"], 
                                     reason_counter["NoDirection"], 
                                     reason_counter["NoAnchor"], 
                                     reason_counter["NoUMI"], 
                                     reason_counter["TrimTooShort"], 
                                     n_pass, 
                                     get_ratio(n_pass, n_total)])) + "\n")
            
    
if __name__ == '__main__':
    main()
    