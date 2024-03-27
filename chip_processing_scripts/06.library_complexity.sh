#!/bin/bash

#$ -N id
#$ -V -cwd
#$ -pe smp 12

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

out_dir=$1
functions_file=$2

source $functions_file
log=${out_dir}/pipeline.log

bam=${out_dir}/Aligned.filt.srt.nodup.bam

prefix=${bam::-4}
pbc_file=${out_dir}/qc/library_complexity.txt

cmd="samtools sort -@ 12 -n -o ${prefix}.srt.tmp.bam ${bam}"
print_exec "$cmd" $log

echo -e "TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF=Distinct/Total\tPBC1=OnePair/Distinct\tPBC2=OnePair/TwoPair" > $pbc_file

bedtools bamtobed -bedpe -i ${prefix}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' >> ${pbc_file}

cmd="rm ${prefix}.srt.tmp.bam"
print_exec "$cmd" $log
