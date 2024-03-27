source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

set -e 

out_dir=$1
functions_file=$2

fastq1=${out_dir}/detect_adapters/reads.1.fastq
fastq2=${out_dir}/detect_adapters/reads.2.fastq
sam=${out_dir}/Aligned.sam
bam=${out_dir}/Aligned.bam
flagstat_qc=${out_dir}/qc/Aligned.flagstat
ref=/reference/public/ucsc/hg38/hg38.fa.gz
log=${out_dir}/pipeline.log

source $functions_file

# 1. Align
cmd="bwa mem -t 12 $ref $fastq1 $fastq2 | samtools sort -@ 12 -O bam -o $bam -"
print_exec "$cmd" $log
       
# 2. Clean
cmd="rm $fastq1 $fastq2"
print_exec "$cmd" $log

# 3. Flagstat
cmd="samtools flagstat -@ 12 ${bam} > ${flagstat_qc}"
print_exec "$cmd" $log
