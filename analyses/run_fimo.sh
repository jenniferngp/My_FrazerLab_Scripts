#!/bin/bash

#$ -N fimo
#$ -V -cwd
#$ -pe smp 1
#$ -o /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/fimo/logs
#$ -e /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/fimo/logs
#$ -t 1-2346:1
#$ -tc 300

# Runs fimo on one motif at a time. Really fast. 

dir=/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/fimo
file=${dir}/memes.txt
seqs=${dir}/sequences.txt
meme=`tail -n +$SGE_TASK_ID $file | head -1`

cmd="fimo --o ${dir}/fimo_out/${meme} --thresh 0.01 --max-stored-scores 1000000 --no-qvalue --skip-matched-sequence --text ${dir}/jaspar_meme/${meme} $seqs > ${dir}/fimo_out/${meme}.txt"
echo $cmd; eval $cmd