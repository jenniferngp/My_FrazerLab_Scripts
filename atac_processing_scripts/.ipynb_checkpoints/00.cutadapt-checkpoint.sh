#!/bin/bash

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

set -e # exit script when an error is encountered

script_dir=/frazer01/home/jennifer/software/encode-atac-seq-pipeline

fq_dir=$1 # fastq directory
uuid=$2 # uuid
out_dir=$3 # output directory
log=$4 # pipeline log
functions_file=$5
cur_dir=`pwd`

source $functions_file

echo "####### Running 00.cutadapt.sh #######" >> $log
echo fq_dir=$fq_dir >> $log
echo uuid=$uuid >> $log
echo out_dir=$out_dir >> $log
echo log=$log >> $log
echo current_directory=$cur_dir >> $log

# 0. Make a scratch directory
if [ ! -d scratch ]
then 
    echo "Missing scratch directory. Created ${cur_dir}/scratch" >> $log
    print_exec "mkdir scratch" $log
fi

dir=`mktemp -d -p ${cur_dir}/scratch`
echo "Created ${dir} to be the working directory" >> $log

# If working directory does not exit, exit script.
if [ ! -d ${dir} ]; then exit 1; fi

# 1. concatenate fastqs
cmd="rsync ${fq_dir}/${uuid}/*.fastq.gz ${dir}"
print_exec "$cmd" $log

cmd="zcat ${dir}/*R1*.fastq.gz > ${dir}/R1.fastq"
print_exec "$cmd" $log

cmd="zcat ${dir}/*R2*.fastq.gz > ${dir}/R2.fastq"
print_exec "$cmd" $log

# 2. remove original fastqs to reduce memory storage
cmd="rm ${dir}/*.fastq.gz"
print_exec "$cmd" $log

in_files=( ${dir}/R1.fastq ${dir}/R2.fastq )

echo Fastq input.. ${in_files[@]} >> $log

# 2. detect which adapters were used (Ilumina, etc.)
#fqs=($fastqs)
for fq in ${in_files[@]}; do
    prefix=`basename $fq .fastq.gz`
	out=${out_dir}/detect_adapters/${prefix}.adapter.txt
	cmd="python3 ${script_dir}/src/detect_adapter.py $fq > $out"
	print_exec "$cmd" $log
done

# 3. trim adapters from the reads using cutadapt
ADAPTER_FWD=`cat ${out_dir}/detect_adapters/R1*.adapter.txt`
ADAPTER_REV=`cat ${out_dir}/detect_adapters/R2*.adapter.txt`
echo "ADAPTER_FWD" $ADAPTER_FWD
echo "ADAPTER_REV" $ADAPTER_REV
cmd="cutadapt --cores 4 -a $ADAPTER_FWD -A $ADAPTER_REV -o ${out_dir}/detect_adapters/reads.1.fastq -p ${out_dir}/detect_adapters/reads.2.fastq ${in_files[0]} ${in_files[1]}"
print_exec "$cmd" $log

# 4. clean
cmd="rm -r $dir"
print_exec "$cmd" $log


