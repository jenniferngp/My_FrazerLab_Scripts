#!/bin/bash

#$ -N stats
#$ -V -cwd
#$ -pe smp 4
#$ -o logs/stats
#$ -e logs/stats

pipe_dir=$1
sample_list=$2
uuid=`tail -n +2 $sample_list | tail -n +$SGE_TASK_ID | cut -f2 | head -1`

source activate encode-chip
set -e

bam=${pipe_dir}/${uuid}/Aligned.merged.bam

# 8. Index bam
date >& 2
if [ ! ${bam}.bai ]; then rm ${bam}.bai; fi
cmd="samtools index -@ 4 ${bam} ${bam}.bai"
echo $cmd >& 2; eval $cmd

# 9. Run flagstat
date >& 2
cmd="samtools flagstat -@ 4 ${bam} > ${bam}.flagstat"
echo $cmd >& 2; eval $cmd

# 10. Run idxstats
date >& 2
cmd="samtools idxstats -@ 4 ${bam} > ${bam}.idxstats"
echo $cmd >& 2; eval $cmd

# 10. Run samtools stats
date >& 2
cmd="samtools stats -@ 4 ${bam} > ${bam}.stats"
echo $cmd >& 2; eval $cmd

# 11. bigWig
cmd="bamCoverage --numberOfProcessors 4 --bam ${bam} --outFileName ${pipe_dir}/${uuid}/Aligned.merged.bw"
echo $cmd >& 2; eval $cmd

date >& 2