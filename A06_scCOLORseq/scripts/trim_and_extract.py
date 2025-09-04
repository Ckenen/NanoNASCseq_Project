#!/usr/bin/env python
import sys
from collections import defaultdict
import pysam
import edlib
import numpy as np

mapper = {"A": "T", "C": "G", "G": "C", "T": "A"}


def reverse_complement(seq):
    return "".join([mapper[b] for b in seq[::-1]])


def align(query, reference):
    a = edlib.align(query, reference, task="locations", mode="HW")
    ed = a["editDistance"]
    start, end = a["locations"][0]
    end += 1
    return start, end, ed


def find_polyt(seq):
    m = np.zeros(len(seq) + 1)
    for i, base in enumerate(seq):
        if base == "T":
            m[i + 1] = m[i] + 1
        else:
            m[i + 1] = max(m[i] - 2, 0)
    vmax = max(m)
    start = 0
    end = None
    for i, v in enumerate(m):
        if v == 0:
            start = i
        if v == vmax:
            end = i
            break
    return start, end


def main():
    infile, prefix = sys.argv[1:]
    
    counter = defaultdict(int)
    unambiguous = 0
    ambiguous = 0
        
    adapter_head = "AAGCAGTGGTATCAACGCAGAGT"
    adapter_tail = "ACTCTGCGTTGATACCACTGCTT"
    
    fw1 = open(prefix + "_unambiguous_barcode_R1.fastq", "w+")
    fw2 = open(prefix + "_unambiguous_barcode_R2.fastq", "w+")
    fw3 = open(prefix + "_ambiguous_barcode_R1.fastq", "w+")
    fw4 = open(prefix + "_ambiguous_barcode_R2.fastq", "w+")
    fw5 = open(prefix + "_whitelist.txt", "w+")
    
    barcodes = set()

    with pysam.FastxFile(infile) as f:
        for i, read in enumerate(f):
            total += 1
            
            seq = read.sequence
            qua = read.quality
            
            # raw length
            
            if len(seq) < 400:
                counter["TooShort"] += 1
                continue
                
            # adapter
            
            w = 100
            seq_head = seq[:w]
            seq_tail = seq[-w:]
            offset = len(seq) - w
            x1, y1, ed1 = align(adapter_head, seq_head)
            x2, y2, ed2 = align(adapter_tail, seq_tail)
            x2, y2 = x2 + offset, y2 + offset
            if ed1 > 6 or ed2 > 6:
                counter["NoAdater"] += 1
                continue
            else:
                seq = seq[y1:x2]
                qua = qua[y1:x2]
                
            # direction
            
            w = 200
            seq_head = seq[:w]
            seq_tail = seq[-w:]
            
            direction = None
            x1, y1 = find_polyt(seq_head)
            x2, y2 = find_polyt("".join([mapper[b] for b in seq_tail[::-1]]))
            count1 = y1 - x1
            count2 = y2 - x2
            
            flag1 = count1 >= 15 and x1 >= 37 and x1 <= 47 # is polyT
            flag2 = count2 >= 15 and x2 >= 37 and x2 <= 47 # is polyA
            if flag1 and not flag2:
                x, y = x1, y1
            elif not flag1 and flag2:
                x, y = x2, y2
                seq = "".join([mapper[b] for b in seq[::-1]])
                qua = qua[::-1]
            else:
                counter["NoDirection"] += 1
                continue
                
            bc_seq = seq[2:26]
            bc_qua = qua[2:26]
            umi_seq = seq[26:42]
            umi_qua = qua[26:42]
            
            rna_seq = seq[y2:]
            rna_qua = qua[y2:]
            rna_seq = "".join([mapper[b] for b in rna_seq[::-1]])
            rna_qua = rna_qua[::-1]

            is_perfect_bc = True
            for j in range(0, len(bc_seq), 2):
                if bc_seq[j] != bc_seq[j + 1]:
                    is_perfect_bc = False
                    break
            
            bc_umi_seq = bc_seq + umi_seq
            bc_umi_qua = bc_qua + umi_qua
            
            if is_perfect_bc:
                barcodes.add(bc_seq)
                fw1.write("@%s\n%s\n+\n%s\n" % (read.name, bc_umi_seq, bc_umi_qua))
                fw2.write("@%s\n%s\n+\n%s\n" % (read.name, rna_seq, rna_qua))
                unambiguous += 1
            else:
                fw3.write("@%s\n%s\n+\n%s\n" % (read.name, bc_umi_seq, bc_umi_qua))
                fw4.write("@%s\n%s\n+\n%s\n" % (read.name, rna_seq, rna_qua))                        
                ambiguous += 1
            counter["Pass"] += 1

    for bc in sorted(barcodes):
        fw5.write("%s\n" % bc)
        
    fw1.close()
    fw2.close()
    fw3.close()
    fw4.close()
    fw5.close() 
    
    print("Total reads:", total)
    print("Too short reads:", counter["TooShort"])
    print("No adapter reads:", counter["NoAdapter"])
    print("No direction reads:", counter["NoDirection"])
    print("Pass reads:", counter["Pass"])
    
    print("Unambiguous reads:", unambiguous)
    print("Ambiguous reads:", ambiguous)
    print("Whitelist barcodes:", len(barcodes))
    

if __name__ == "__main__":
    main()