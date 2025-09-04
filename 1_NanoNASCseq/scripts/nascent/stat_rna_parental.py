#!/usr/bin/env python
import sys
from collections import Counter, defaultdict
import pysam


def main():
    bamfile, txtfile = sys.argv[1:]
    
    parentals = defaultdict(list)
    with pysam.AlignmentFile(bamfile) as f:
        for s in f:
            cn = s.get_tag("CN")
            pm = s.get_tag("PM")
            parentals[cn].append(pm)
            
    parental = dict()
    for k, v in parentals.items():
        counter = Counter(v)
        items = list(sorted(counter.items(), key=lambda item: item[1]))
        pm = items[-1][0]
        if len(items) >= 2 and items[-1][1] < items[-2][1] * 2:
            pm = "U"
        parental[k] = pm
    
    with open(txtfile, "w+") as fw:
        for k, v in parental.items():
            fw.write("%s\t%s\n" % (k, v))
            
            
if __name__ == "__main__":
    main()
        
        