#!/usr/bin/env bash
gtf1=$1 # query.gtf
gtf2=$2 # reference.gtf
fasta=$3 # reference.fasta
threads=$4 # threads
outdir=$5 # output directory
name=`basename ${gtf1} .gtf`

export PYTHONPATH="${HOME}/software/cDNA_Cupcake-27.0.0:${PYTHONPATH}"
export PYTHONPATH="${HOME}/software/cDNA_Cupcake-27.0.0/sequence:${PYTHONPATH}"
export PATH="${HOME}/software/SQANTI3-4.2:${PATH}"

sqanti3_qc.py -d ${outdir} --cpus ${threads} --report pdf ${gtf1} ${gtf2} ${fasta} && \
sqanti3_RulesFilter.py --report skip ${outdir}/${name}_classification.txt ${outdir}/${name}_corrected.fasta ${outdir}/${name}_corrected.gtf