#!/bin/sh
set -e

bam=$1
prefix=$2

# Outputs:
# ${prefix}_raw_both.bw
# ${prefix}_raw_pos.bw
# ${prefix}_raw_neg.bw
# ${prefix}_norm_both.bw
# ${prefix}_norm_pos.bw
# ${prefix}_norm_neg.bw

bedtools bamtobed -bed12 -split -i ${bam} > ${prefix}.bed
wc -l ${prefix}.bed | awk '{print $1}' > ${prefix}_lib_size.txt
samtools view -H ${bam} | grep '@SQ' | awk -v OFS='\t' '{print $2,$3}' | sed 's/SN://g' | sed 's/LN://g' > ${prefix}_genome_size.txt

lib_size=`cat ${prefix}_lib_size.txt`
scale=`echo ${lib_size} | awk '{print 1000000/$1}'`

# raw

bedtools genomecov -bg -split -i ${prefix}.bed -g ${prefix}_genome_size.txt | sort -k1,1 -k2,2n > ${prefix}_raw_both.bedGraph
bedGraphToBigWig ${prefix}_raw_both.bedGraph ${prefix}_genome_size.txt ${prefix}_raw_both.bw

bedtools genomecov -bg -split -strand + -i ${prefix}.bed -g ${prefix}_genome_size.txt | sort -k1,1 -k2,2n > ${prefix}_raw_pos.bedGraph
bedGraphToBigWig ${prefix}_raw_pos.bedGraph ${prefix}_genome_size.txt ${prefix}_raw_pos.bw

bedtools genomecov -bg -split -strand - -i ${prefix}.bed -g ${prefix}_genome_size.txt | sort -k1,1 -k2,2n > ${prefix}_raw_neg.bedGraph
bedGraphToBigWig ${prefix}_raw_neg.bedGraph ${prefix}_genome_size.txt ${prefix}_raw_neg.bw

# normalized

bedtools genomecov -bg -split -scale ${scale} -i ${prefix}.bed -g ${prefix}_genome_size.txt | sort -k1,1 -k2,2n > ${prefix}_norm_both.bedGraph
bedGraphToBigWig ${prefix}_norm_both.bedGraph ${prefix}_genome_size.txt ${prefix}_norm_both.bw

bedtools genomecov -bg -split -scale ${scale} -strand + -i ${prefix}.bed -g ${prefix}_genome_size.txt | sort -k1,1 -k2,2n > ${prefix}_norm_pos.bedGraph
bedGraphToBigWig ${prefix}_norm_pos.bedGraph ${prefix}_genome_size.txt ${prefix}_norm_pos.bw

bedtools genomecov -bg -split -scale ${scale} -strand - -i ${prefix}.bed -g ${prefix}_genome_size.txt | sort -k1,1 -k2,2n > ${prefix}_norm_neg.bedGraph
bedGraphToBigWig ${prefix}_norm_neg.bedGraph ${prefix}_genome_size.txt ${prefix}_norm_neg.bw

rm ${prefix}.bed
rm ${prefix}_lib_size.txt
rm ${prefix}_genome_size.txt
rm ${prefix}_raw_both.bedGraph
rm ${prefix}_raw_pos.bedGraph
rm ${prefix}_raw_neg.bedGraph
rm ${prefix}_norm_both.bedGraph
rm ${prefix}_norm_pos.bedGraph
rm ${prefix}_norm_neg.bedGraph

