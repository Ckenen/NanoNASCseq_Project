#!/bin/sh
fq=$1
echo -e "File\tReads\tBases"
zcat ${fq} | awk -v OFS='\t' '{if(NR%4==2){a+=1;b+=length($1)}}END{print "'${fq}'",a,b}'
