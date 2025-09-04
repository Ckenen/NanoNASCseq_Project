#!/usr/bin/env bash
fastq=$1
echo -e 'Length\tCount'
gzip -d -c $fastq | awk '{if(NR%4==2){print length($1)}}' | sort -k1,1 | uniq -c | awk '{print $2"\t"$1}' | sort -k1,1n