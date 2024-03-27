#!/bin/bash

#$ -N qc
#$ -V -cwd
#$ -pe smp 8

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

set -e

# ATAC-seq pipeline
# Documentation: https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
# Github: https://github.com/ENCODE-DCC/atac-seq-pipeline/tree/master

# ===================
# Create tagAlign file
# ===================
out_dir=$1
log=$2
functions_file=$3

prefix=${out_dir}/Aligned.sorted
final_bam_prefix=${out_dir}/Aligned.sorted.filt.nodup.nomito
final_bam_file=${final_bam_prefix}.bam
final_bedpe_file=${final_bam_prefix}.bedpe.gz
final_ta_file=${final_bam_file}.tagAlign.gz
script_dir=/projects/PPC/pipeline/ATAC-Seq/work/2023_0720/scripts
pbc_file_qc=${final_bam_prefix}.pbc.qc

source $functions_file

echo "####### Running 06.qc.sh #######" >> $log

# 1. sort bam file by name
cmd="samtools sort -@ 8 -n ${final_bam_file} -o ${final_bam_prefix}.srt.tmp.bam"
print_exec "$cmd" $log

# 2. compute library complexity

# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

bedtools bamtobed -bedpe -i ${final_bam_prefix}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${pbc_file_qc}

cmd="mv ${pbc_file_qc} ${out_dir}/qc"
print_exec "$cmd" $log

cmd="rm ${final_bam_prefix}.srt.tmp.bam"
print_exec "$cmd" $log

# 3. convert bam to bedpe file 
cmd="bedtools bamtobed -bedpe -mate1 -i ${final_bam_prefix}.srt.tmp.bam | gzip -nc > ${final_bedpe_file}"
print_exec "$cmd" $log

# 4. restructure to tagalign file
zcat ${final_bedpe_file} | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${final_ta_file}

# 5. apply tn5 shifting to tag align file
shifted_tag_file="${final_bam_prefix}.tn5.tagAlign.gz"

zcat -f ${final_ta_file} | awk 'BEGIN {OFS = "\t"} { if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} if ($2 >= $3) { if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1} } print $0}' | gzip -nc > ${shifted_tag_file}

# 6. compute insert size metrics
insert_data=${final_bam_prefix}.insertsize.metrics
insert_plot=${final_bam_prefix}.insertsize.metrics.plot

cmd="picard CollectInsertSizeMetrics \
INPUT=${final_bam_file} OUTPUT=${insert_data} H=${insert_plot} \
VERBOSITY=ERROR QUIET=TRUE \
USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE \
W=1000 STOP_AFTER=5000000"
print_exec "$cmd" $log

cmd="python3 ${script_dir}/plot_fragment_length.py $insert_data $final_bam_prefix"
print_exec "$cmd" $log

cmd="mv $insert_data ${out_dir}/*.plot ${out_dir}/*.png ${out_dir}/*.qc ${out_dir}/qc"
print_exec "$cmd" $log

cmd="rm $final_bedpe_file $final_ta_file $shifted_tag_file"
print_exec "$cmd" $log
