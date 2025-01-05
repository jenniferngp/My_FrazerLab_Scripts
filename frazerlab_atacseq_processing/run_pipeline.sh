#!/bin/bash

prj_dir=/projects/PPC
data_dir=${prj_dir}/data/PPC_Pilot/ATAC-seq
pipe_dir=${prj_dir}/pipeline/ATAC-Seq/sample
script_dir=/projects/PPC/pipeline/ATAC-Seq/work/2021_0909/script

samples_file=$1

id=`tail -n +$SGE_TASK_ID $samples_file | head -1`

if [ -d ${pipe_dir}/${id} ]; then exit; fi

mkdir ${pipe_dir}/${id}

cmd="${script_dir}/star.pbs ${data_dir}/${id} ${pipe_dir}/${id} $id"
qsub -N star_${id} $cmd

cmd="${script_dir}/filt.pbs ${pipe_dir}/${id}/Aligned.out.bam ${pipe_dir}/${id}"
qsub -N filt_${id} -hold_jid star_${id} -l short $cmd

cmd="${script_dir}/identity.pbs ${pipe_dir}/${id}/Aligned.out.filt.bam ${pipe_dir}/${id}"
qsub -N ident_${id} -hold_jid filt_${id} $cmd

cmd="${script_dir}/peak.pbs ${pipe_dir}/${id}/Aligned.out.filt.bam broad ${pipe_dir}/${id}"
qsub -N peak_${id} -hold_jid filt_${id} -l short $cmd

cmd="${script_dir}/peak.pbs ${pipe_dir}/${id}/Aligned.out.filt.bam narrow ${pipe_dir}/${id}"
qsub -N peak_${id} -hold_jid filt_${id} -l short $cmd

cmd="${script_dir}/reads_in_promoter.pbs ${pipe_dir}/${id} ${pipe_dir}/${id}"
qsub -N prom_${id} -hold_jid filt_${id} -l short $cmd

cmd="${script_dir}/count.pbs ${pipe_dir}/${id}/Aligned.out.filt.bam ${pipe_dir}/${id}/broad_peaks.broadPeak ${pipe_dir}/${id}/broad_peaks.count"
qsub -N count_${id} -hold_jid peak_${id} -l short $cmd

cmd="${script_dir}/fsize.pl ${pipe_dir}/${id}/Aligned.out.filt.bam > ${pipe_dir}/${id}/Aligned.out.filt.fs"
echo "date >& 2; echo $cmd >& 2; $cmd; date >& 2" |\
qsub -V -cwd -e logs -o logs -N fsize_${id} -hold_jid filt_${id} -l short

cmd="samtools view -c -q 20 ${pipe_dir}/${id}/Aligned.out.mdup.bam > ${pipe_dir}/${id}/Aligned.out.mdup.q20"
echo "date >& 2; echo '$cmd' >& 2; $cmd; date >& 2" |\
qsub -V -cwd -e logs -o logs -N q20_${id} -hold_jid filt_${id} -l short

