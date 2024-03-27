#!/bin/bash

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

set -e

ENCODE_DIR=/frazer01/home/jennifer/software/encode-atac-seq-pipeline

out_dir=$1
log=$2
functions_file=$3
fastq1=${out_dir}/detect_adapters/reads.1.fastq
fastq2=${out_dir}/detect_adapters/reads.2.fastq
bam_file=${out_dir}/Aligned.bam
bam_sorted_file=${out_dir}/Aligned.sorted.bam
flagstat_file=${out_dir}/qc/Aligned.sorted.flagstat

ref=/reference/public/ucsc/hg38/hg38.fa.gz

source $functions_file

echo "####### Running 01.align.sh #######" >> $log

# old command for bowtie2 alignment
#BWT_IDX=/reference/private/Gencode.v34lift37/bowtie2_index/bowtie2_index
#cmd="bowtie2 -k 4 -X 2000 --mm --threads 8 -x $BWT_IDX -1 $FASTQ1 -2 $FASTQ2 > ${SAM_FILE}"
#echo $cmd; eval $cmd

# 1. align with bwa-mem, which outputs a sam file. then, use samtools to convert to bam
cmd="bwa mem -t 8 $ref $fastq1 $fastq2 | samtools view -@ 12 -O bam -o $bam_file -"
print_exec "$cmd" $log

# 2. sort bam file
cmd="samtools sort -@ 8 -o $bam_sorted_file $bam_file"
print_exec "$cmd" $log

# 3. remove fastq files to reduce storage
cmd="rm $fastq1 $fastq2 $bam_file"
print_exec "$cmd" $log

# 4. generate flagstat metrics
cmd="samtools flagstat -@ 8 $bam_sorted_file > $flagstat_file"
print_exec "$cmd" $log

# 5. index bam file
cmd="samtools index -@ 8 $bam_sorted_file"
print_exec "$cmd" $log





