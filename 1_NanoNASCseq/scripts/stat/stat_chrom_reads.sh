#!/bin/sh

bam=$1
out=$2

samtools view -@ 4 -m 200 -q 30 -F 2308 $bam \
    | awk '{print $3}' | sort | uniq -c \
    | awk -v OFS='\t' '{print $2,$1}' | sort -k1,1 > $out