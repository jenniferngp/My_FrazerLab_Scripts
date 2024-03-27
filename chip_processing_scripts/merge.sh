#!/bin/bash

#$ -N merge
#$ -o logs/merge
#$ -e logs/merge
#$ -V -cwd
#$ -pe smp 8

merge_map=$1
pipe_dir=$2
line=(`tail -n +2 $merge_map | tail -n +$SGE_TASK_ID | head -1`)
merge_uuid=${line[3]}
npremerge=${line[4]}
out_dir=${pipe_dir}/${merge_uuid}
log_file=${out_dir}/merge.log

echo "Merge uuid:" $merge_uuid
echo "Number of premerges:" $npremerge
echo "Output directory:" $out_dir
echo "Log file:" $log_file

if [ ! -d $out_dir ]; then mkdir $out_dir; fi

# 1. merge
if [ $npremerge == 1 ]
then
    cmd="samtools merge -o ${out_dir}/Aligned.merged.bam ${pipe_dir}/${line[0]}/Aligned.filt.srt.nodup.bam"
elif [ $npremerge == 2 ]
then
    cmd="samtools merge -o ${out_dir}/Aligned.merged.bam ${pipe_dir}/${line[0]}/Aligned.filt.srt.nodup.bam ${pipe_dir}/${line[1]}/Aligned.filt.srt.nodup.bam"
elif [ $npremerge == 3 ]
then
    cmd="samtools merge -o ${out_dir}/Aligned.merged.bam ${pipe_dir}/${line[0]}/Aligned.filt.srt.nodup.bam ${pipe_dir}/${line[1]}/Aligned.filt.srt.nodup.bam ${pipe_dir}/${line[2]}/Aligned.filt.srt.nodup.bam"
fi

date >> $log_file

if [ -f ${out_dir}/Aligned.merged.bam ]
then 
    echo "Already done!" >> $log_file
    echo "Already done!" >& 2
    exit 1
fi

echo $cmd >> $log_file

echo $cmd >& 2; eval $cmd
