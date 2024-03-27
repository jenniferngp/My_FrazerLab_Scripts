#!/bin/bash

#$ -N peaks
#$ -V -cwd
#$ -o logs/peaks
#$ -e logs/peaks
#$ -pe smp 12

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

set -e

bam=$1
name=${bam::-4}
functions_file=$2

source $functions_file

# 1. sort bam file
cmd="samtools sort -@ 12 -n -o ${name}.nsort.bam ${bam}"
echo $cmd >& 2
eval $cmd

# 2. fix mate information
cmd="samtools fixmate -@ 12 ${name}.nsort.bam ${name}.nsort.fixed.bam"
echo $cmd >& 2
eval $cmd

# 3. convert bam to bed file
cmd="bedtools bamtobed -bedpe -mate1 -i ${name}.nsort.fixed.bam | gzip -nc > ${name}.tmp"
echo $cmd >& 2
eval $cmd

# 4. create tag align file, which is a bed file
zcat ${name}.tmp | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${name}.tagAlign.gz

cmd="rm ${name}.nsort.bam ${name}.nsort.fixed.bam ${name}.tmp"
echo $cmd >& 2
eval $cmd

echo "Saved: ${name}.tagAlign.gz" >& 2
