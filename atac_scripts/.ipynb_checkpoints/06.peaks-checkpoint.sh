#!/bin/bash

#$ -N peaks
#$ -V -cwd
#$ -pe smp 4

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

set -e

out_dir=$1
log=$2
functions_file=$3

prefix=${out_dir}/peaks/narrow
pval_thresh=0.01 # default in ENCODE
NPEAKS=300000 # capping number of peaks called from MACS2c
#smooth_window=150 # default in ENCODE
#shiftsize=$(( -$smooth_window/2 ))
#gensz=hs # human genome size
smooth_window=$4 #150
shiftsize=$5 #75
gensz=$6 
script_dir=/projects/PPC/pipeline/ATAC-Seq/work/2023_0720/scripts

source $functions_file

echo "####### Running 07.peaks.sh #######" >> $log

if [ ! -d ${out_dir}/peaks ]; then mkdir ${out_dir}/peaks; fi

# 1. call narrow peaks using macs2
bam=${out_dir}/Aligned.sorted.filt.nodup.nomito.bam
cmd="macs2 callpeak -t $bam -f BAMPE -n ${out_dir}/peaks/narrow -g $gensz -q $pval_thresh --shift $shiftsize  --extsize $smooth_window --nomodel -B --SPMR --keep-dup all --call-summits"
print_exec "$cmd" $log

# 2. remove peaks that overlap blacklist regions
blacklist=/reference/public/ENCODE/hg38-blacklist.v2.bed
filtered_peak=${out_dir}/peaks/narrow_peaks_noblacklist.narrowPeak
bedtools intersect -v -a ${out_dir}/peaks/narrow_peaks.narrowPeak -b ${blacklist} | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' | grep -P 'chr[\dXY]+[ \t]'   > ${filtered_peak}

# 3. run feature counts to get the number of reads per peak. Keep the log to obtain the FRIP metric.
bam=${out_dir}/Aligned.sorted.filt.nodup.nomito.bam
peak=${out_dir}/peaks/narrow_peaks_noblacklist.narrowPeak

# convert bed to saf
cmd="cat $peak | ${script_dir}/peaks2saf.pl > ${out_dir}/saf"
print_exec "$cmd" $log

# run feature counts
cmd="featureCounts -p --countReadPairs -B -C -T 4 -F SAF -a ${out_dir}/saf -o ${peak}.counts ${bam}"
print_exec "$cmd" $log

# clean
cmd="rm ${out_dir}/saf"
print_exec "$cmd" $log
