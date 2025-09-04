#!/usr/bin/env python
import sys
from Bio.Seq import Seq

for i, line in enumerate(sys.stdin):
    j = i % 4
    if j == 1:
        line = str(Seq(line[:-1]).reverse_complement()) + "\n"
    elif j == 3:
        line = line[:-1][::-1] + "\n"
    sys.stdout.write(line)