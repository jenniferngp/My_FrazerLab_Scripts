#!/bin/bash

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

set -e

out_dir=$1
log=$2
functions_file=$3
bam_file=${out_dir}/Aligned.sorted.bam
mito_flagstat=${out_dir}/qc/Aligned.sorted.mito.flagstat
mito_bam=${out_dir}/Aligned.sorted.mito.bam

source $functions_file

echo "####### Running 02.mito.sh #######" >> $log

# 1. create bam file that only contains mitochondrial reads
cmd="samtools view -b ${bam_file} chrM > ${mito_bam}"
print_exec "$cmd" $log

# 2. run flagstat on mito bam
cmd="samtools flagstat -@ 8 ${mito_bam} > ${mito_flagstat}"
print_exec "$cmd" $log

# 3. remove mito bam
cmd="rm $mito_bam"
print_exec "$cmd" $log

