#!/bin/bash

#$ -V -cwd
#$ -pe smp 8

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

set -e

encode_dir=/frazer01/home/jennifer/software/encode-atac-seq-pipeline

uuid=$1
fq_dir=$2
out_dir=$3
functions_file=$4

echo $functions_file

source $functions_file

log=${out_dir}/pipeline.log

fastqs="${out_dir}/detect_adapters/R1.fastq.gz ${out_dir}/detect_adapters/R2.fastq.gz"

# 0. Concatenate R1 fastq from all sequencing runs
cmd="cat ${fq_dir}/*R1* > ${out_dir}/detect_adapters/R1.fastq.gz"
print_exec "$cmd" $log

# 0. Concatenate R2 fastq from all sequencing runs
cmd="cat ${fq_dir}/*R2* > ${out_dir}/detect_adapters/R2.fastq.gz"
print_exec "$cmd" $log

# 1. Detect which adapters were used 
fqs=($fastqs)
for fq in ${fqs[@]}; do
    prefix=`basename $fq .fastq.gz`
    out=${out_dir}/detect_adapters/${prefix}.adapter.txt
    cmd="python3 ${encode_dir}/src/detect_adapter.py $fq > $out"
    print_exec "$cmd" $log
done

# 2. Trim adapters from reads using cutadapt
ADAPTER_FWD=`cat ${out_dir}/detect_adapters/*R1*.adapter.txt`
ADAPTER_REV=`cat ${out_dir}/detect_adapters/*R2*.adapter.txt`

echo "ADAPTER_FWD:" $ADAPTER_FWD >& 2
echo "ADAPTER_REV:" $ADAPTER_REV >& 2

cmd="cutadapt --cores 12 -a $ADAPTER_FWD -A $ADAPTER_REV -o ${out_dir}/detect_adapters/reads.1.fastq -p ${out_dir}/detect_adapters/reads.2.fastq ${fqs[0]} ${fqs[1]}"
print_exec "$cmd" $log

# 3. Run fastqc
cmd="fastqc -o ${out_dir}/detect_adapters ${out_dir}/detect_adapters/reads.*.fastq"
print_exec "$cmd" $log

# 4. Unzip fastqc
cmd="unzip ${out_dir}/detect_adapters/reads.1_fastqc.zip -d ${out_dir}/detect_adapters; unzip ${out_dir}/detect_adapters/reads.2_fastqc.zip -d ${out_dir}/detect_adapters"

# 5. Clean to clear space
cmd="rm ${out_dir}/detect_adapters/R1.fastq.gz ${out_dir}/detect_adapters/R2.fastq.gz"
print_exec "$cmd" $log