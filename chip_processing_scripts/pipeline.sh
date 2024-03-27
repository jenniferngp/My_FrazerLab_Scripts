#!/bin/bash/

#$ -N chip
#$ -t 292-292:1
#$ -tc 5
#$ -pe smp 12
#$ -V -cwd
#$ -o logs
#$ -e logs

# qsub pipeline [uuid_list] [fq_dir] [out_dir] [functions_file] [script_dir] [upperbound] [call_peaks] [input]

uuid_list=$1
fq_dir=$2
out_dir=$3
functions_file=$4
script_dir=$5
upperbound=$6 # 50 for TF, 100 for histone
call_peaks=$7
input=$8

uuid=`tail -n +$SGE_TASK_ID $uuid_list | head -1`
out_dir=${out_dir}/${uuid}

source $functions_file
source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

# get genome reference
# (already done, do not re-run)
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.
# bwa=/software/bwa-0.7.17/bwa
# ${bwa} index hg38.fa.gz       

mkdir ${out_dir} ${out_dir}/detect_adapters ${out_dir}/qc

# 0. Concatenate fastqs. Trim adapters. Run QC on fastq files. 
sh ${script_dir}/00.pre-align.sh $uuid ${fq_dir}/${uuid} $out_dir $functions_file

# 1. Align
sh ${script_dir}/01.align.sh $out_dir $functions_file

# 2. Filter reads
sh ${script_dir}/02.filter.sh $out_dir $functions_file

# 3. Create tagAlign
sh ${script_dir}/03.tagAlign.sh ${out_dir}/Aligned.filt.srt.nodup.bam $functions_file

# 4. Cross correlation
sh ${script_dir}/04.cross_correlation.sh $out_dir $upperbound $functions_file

# 5. Identify
sh ${script_dir}/05.identity.sh $out_dir $functions_file

# 6. Library complexity
sh ${script_dir}/06.library_complexity.sh $out_dir $functions_file

# 7. Call peaks
if [ $call_peaks == T ]
then
    input_dir=${out_dir}/${input}
    
    if [ ! -d ${out_dir}/logs ]; then mkdir ${out_dir}/logs; fi
    
    log_out=${out_dir}/logs/peaks.out
    log_err=${out_dir}/logs/peaks.err
    
    qsub -N peaks -V -cwd -pe smp 4 -o $log_out -e $log_err ${script_dir}/07.peaks.sh ${out_dir} ${input_dir} ${functions_file}
fi

# 8. Clean
rm ${out_dir}/Aligned.bam 
