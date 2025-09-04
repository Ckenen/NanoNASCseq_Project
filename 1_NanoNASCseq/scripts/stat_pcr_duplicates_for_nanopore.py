#!/usr/bin/env python
import sys
import edlib
import pandas as pd
from pyBioInfo.IO.File import BamFile
from pyBioInfo.Utils import BundleBuilder


def cmp_seq(seq, ref):
    return edlib.align(seq, ref, task='path', mode='HW')


def main():
    infile, outfile = sys.argv[1:]

    alignments = []
    with BamFile(infile) as f:
        for align in f:
            assert align.name[-9] == ":"
            align.umi = align.name[-8:]
            alignments.append(align)
            
    bundles = []
    for bundle in BundleBuilder(alignments, keep=True):
        bundles.append(bundle)
        
    rows = []
    for max_polya_diff in range(0, 10):
        for max_edit_distance in range(0, 5):
            total_groups = []
            print(max_polya_diff, max_edit_distance, sep="\t")
            for bundle in bundles:
                array1 = [] # +
                array2 = [] # -
                for obj in bundle.data:
                    if obj.strand == "+":
                        array1.append(obj)
                    else:
                        array2.append(obj)
                array1 = list(sorted(array1, key=lambda item: item.end))
                array2 = list(sorted(array2, key=lambda item: item.start))

                groups = []

                # +
                tmp = array1.copy()
                while len(tmp) > 0:
                    obj1 = tmp[0]
                    tmp.pop(0)
                    group = [obj1]
                    j = 0
                    while j < len(tmp):
                        obj2 = tmp[j]
                        if obj2.end - obj1.end <= max_polya_diff:
                            a = cmp_seq(obj1.umi, obj2.umi)
                            if a["editDistance"] <= max_edit_distance:
                                group.append(obj2)
                                tmp.pop(j)
                            else:
                                j += 1
                        else:
                            break
                    groups.append(group)

                # -
                tmp = array2.copy()
                while len(tmp) > 0:
                    obj1 = tmp[0]
                    tmp.pop(0)
                    group = [obj1]
                    j = 0
                    while j < len(tmp):
                        obj2 = tmp[j]
                        if obj2.start - obj1.start <= max_polya_diff:
                            a = cmp_seq(obj1.umi, obj2.umi)
                            if a["editDistance"] <= max_edit_distance:
                                group.append(obj2)
                                tmp.pop(j)
                            else:
                                j += 1
                        else:
                            break
                    groups.append(group)

                for group in groups:
                    total_groups.append(groups)
                    
            row = [max_polya_diff, max_edit_distance, len(alignments), len(total_groups)]
            rows.append(row)
    
    dat = pd.DataFrame(rows)
    dat.columns = ["MaxPolyADiff", "MaxEditDistance", "Alignments", "Groups"]
    dat["Proportion"] = dat["Groups"] / dat["Alignments"]
    dat.to_csv(outfile, sep="\t", index=False)
    
    
if __name__ == '__main__':
    main()
    