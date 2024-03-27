#!/bin/bash

#$ -pe smp 8
#$ -V -cwd

set -e 

out_dir=$1
log=$2

bam=${out_dir}/Aligned.sorted.filt.nodup.bam
nomito=${out_dir}/Aligned.sorted.filt.nodup.nomito.bam
mito=${out_dir}/Aligned.sorted.filt.nohup.mito.bam

# 1. Make mito-free bam
date >& 2

cmd="samtools idxstats $bam | cut -f 1 | grep -v -P "^chrM" | xargs samtools view $bam -@ 8 -b > ${nomito}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2
    
cmd="samtools index -@ 8 ${nomito} ${nomito}.bai"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2

# 2. Flagstat on mito-free bam
cmd="samtools flagstat -@ 8 ${nomito} > ${out_dir}/qc/Aligned.sorted.filt.nodup.nomito.flagstat.qc"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2

# 3. Make mito bam
cmd="samtools view -b ${bam} chrM > ${mito}"
#echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="samtools flagstat -@ 8 ${mito} > ${out_dir}/qc/Aligned.sorted.filt.nodup.mito.flagstat.qc"
#echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="rm $mito"
#echo $cmd >& 2; echo $cmd >> $log; eval $cmd