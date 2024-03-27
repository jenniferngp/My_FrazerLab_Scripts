#!/bin/bash

#$ -N peaks
#$ -V -cwd
#$ -o logs/peaks
#$ -e logs/peaks
#$ -pe smp 4

# macs2 parameter description: https://manpages.ubuntu.com/manpages/xenial/man1/macs2_callpeak.1.html

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

set -e

out_dir=$1
input_dir=$2
functions_file=$3

log=${out_dir}/macs2.log

source $functions_file

# 0. Get fragment length
cc_qc_log=${out_dir}/qc/Aligned.merged.subsampled.tagAlign.cc.qc
fraglen=`tail -n +2 $cc_qc_log |cut -f3`

# 0. Setup input
bam_chip=${out_dir}/Aligned.filt.srt.nodup.bam
bed_chip=${out_dir}/Aligned.filt.srt.nodup.tagAlign

bam_input=${input_dir}/Aligned.filt.srt.nodup.bam
bed_input=${input_dir}/Aligned.filt.srt.nodup.tagAlign

blacklist=/reference/public/ENCODE/hg38-blacklist.v2.bed

if [ ! -d ${out_dir}/peaks ]; then mkdir ${out_dir}/peaks; fi

# 1. Narrow peaks (bed)
prefix=${out_dir}/peaks/narrow_bed
cmd="macs2 callpeak -t ${bed_chip}.gz -c ${bed_input}.gz -g 3.0e9 -n $prefix -q 1e-2 --nomodel --shift 0 --extsize $fraglen --keep-dup all --call-summits"
print_exec "$cmd" $log

# Remove blacklist peaks
filtered_peak=${prefix}_peaks_noblacklist.narrowPeak
bedtools intersect -v -a ${prefix}_peaks.narrowPeak -b ${blacklist} | grep -P 'chr[\dXY]+[ \t]'   > ${filtered_peak}

# Run featureCounts
cmd="awk -F \"\t\" '{print \$4,\$1,\$2,\$3,\$6}' ${prefix}_peaks_noblacklist.narrowPeak | tr ' ' '\\t' > ${out_dir}/narrow_saf"
print_exec "$cmd" $log

cmd="featureCounts -p --countReadPairs -B -C -T 4 -F SAF -a ${out_dir}/narrow_saf -o ${prefix}_peaks_noblacklist.narrowPeak.counts $bam_chip"
print_exec "$cmd" $log

# 2. Broad peaks (bed)
prefix=${out_dir}/peaks/broad_bed
cmd="macs2 callpeak -t ${bed_chip}.gz -c ${bed_input}.gz -g 3.0e9 -n $prefix -q 1e-2 --nomodel --shift 0 --extsize $fraglen --keep-dup all --broad"
print_exec "$cmd" $log

# Remove blacklist peaks
filtered_peak=${prefix}_peaks_noblacklist.broadPeak
bedtools intersect -v -a ${prefix}_peaks.broadPeak -b ${blacklist} | grep -P 'chr[\dXY]+[ \t]'   > ${filtered_peak}

# Run featureCounts
cmd="awk -F \"\t\" '{print \$4,\$1,\$2,\$3,\$6}' ${prefix}_peaks_noblacklist.broadPeak | tr ' ' '\\t' > ${out_dir}/broad_saf"
print_exec "$cmd" $log

cmd="featureCounts -p --countReadPairs -B -C -T 4 -F SAF -a ${out_dir}/broad_saf -o ${prefix}_peaks_noblacklist.broadPeak.counts $bam_chip"
print_exec "$cmd" $log

# 3. Narrow peaks (bam)
prefix=${out_dir}/peaks/narrow_bam
cmd="macs2 callpeak -t ${bam_chip} -c ${bam_input} -g 3.0e9 -n $prefix -q 1e-2 --nomodel --shift 0 --extsize $fraglen --keep-dup all --call-summits"
print_exec "$cmd" $log

# Remove blacklist peaks
filtered_peak=${prefix}_peaks_noblacklist.narrowPeak
bedtools intersect -v -a ${prefix}_peaks.narrowPeak -b ${blacklist} | grep -P 'chr[\dXY]+[ \t]'   > ${filtered_peak}

# Run featureCounts
cmd="awk -F \"\t\" '{print \$4,\$1,\$2,\$3,\$6}' ${prefix}_peaks_noblacklist.narrowPeak | tr ' ' '\\t' > ${out_dir}/narrow_saf"
print_exec "$cmd" $log

cmd="featureCounts -p --countReadPairs -B -C -T 4 -F SAF -a ${out_dir}/narrow_saf -o ${prefix}_peaks_noblacklist.narrowPeak.counts $bam_chip"
print_exec "$cmd" $log

# 4. Broad peaks (bam)
prefix=${out_dir}/peaks/broad_bam
cmd="macs2 callpeak -t ${bam_chip} -c ${bam_input} -g 3.0e9 -n $prefix -q 1e-2 --nomodel --shift 0 --extsize $fraglen --keep-dup all --broad"
print_exec "$cmd" $log

# Remove blacklist peaks
filtered_peak=${prefix}_peaks_noblacklist.broadPeak
bedtools intersect -v -a ${prefix}_peaks.broadPeak -b ${blacklist} | grep -P 'chr[\dXY]+[ \t]'   > ${filtered_peak}

# Run featureCounts
cmd="awk -F \"\t\" '{print \$4,\$1,\$2,\$3,\$6}' ${prefix}_peaks_noblacklist.broadPeak | tr ' ' '\\t' > ${out_dir}/broad_saf"
print_exec "$cmd" $log

cmd="featureCounts -p --countReadPairs -B -C -T 4 -F SAF -a ${out_dir}/broad_saf -o ${prefix}_peaks_noblacklist.broadPeak.counts $bam_chip"
print_exec "$cmd" $log
