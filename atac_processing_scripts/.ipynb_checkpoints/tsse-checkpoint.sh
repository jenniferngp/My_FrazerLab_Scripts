#!/bin/bash

#$ -N peaks
#$ -V -cwd
#$ -pe smp 4

date >& 2 

set -e

out_dir=$1
log=$2

script_dir=/projects/PPC/pipeline/ATAC-Seq/work/2023_0720/scripts

# TSS enrichment
source activate encode-chip
source activate encode-atac
source activate py2 # for some reason, this only works after you activate encode-chip and encode-atac consecutively

chrsz=/reference/public/hg38/hg38.size.txt
tss_bed=/reference/private/Gencode.v44lift38/promoters.bed 
bam=${out_dir}/Aligned.sorted.filt.nodup.bam
out_prefix=${out_dir}/peaks/Aligned.sorted.filt.nodup
read_len=`cat ${out_dir}/read_length.txt`

cmd="python ${script_dir}/tss_enrichment.py $bam $tss_bed $out_prefix $chrsz $read_len"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd