#!/usr/bin/env python
import sys
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder

f_gtf, f_bed = sys.argv[1:]

with GtfFile(f_gtf) as f, open(f_bed, "w+") as fw:
    records = [x for x in f]
    objs = GtfTranscriptBuilder(records)
    for obj in sorted(objs):
        fw.write(obj.format("BED") + "\n")