#!/usr/bin/env python
import sys
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import BamFile
from nasctools.utils import decode_blocks_str

inbam = sys.argv[1]

with BamFile(inbam, random=False) as f:
    for align in f:
        segment = align.segment
        if segment.is_duplicate:
            continue
        blocks = decode_blocks_str(segment.get_tag("XE"))
        name = segment.get_tag("CN")
        obj = GRange(chrom=align.chrom, blocks=blocks, name=name, strand=align.strand)
        obj.score = segment.get_tag("CS")
        print(obj.format("BED"))
