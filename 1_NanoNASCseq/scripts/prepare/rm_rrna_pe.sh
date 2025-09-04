#!/bin/sh
set -e

fq1=$1
fq2=$2
bt2idx=$3
threads=$4
prefix=$5

bowtie2 -p ${threads} --local --no-unal --un-conc ${prefix}.fq -x ${bt2idx}/ref -1 ${fq1} -2 ${fq2} \
    | samtools view -@ ${threads} -b - \
    | samtools sort -@ ${threads} -T ${prefix}_TMP - > ${prefix}.bam
pigz -p ${threads} ${prefix}.1.fq ${prefix}.2.fq
