#!/bin/bash

prj_dir=/projects/PPC
data_dir=${prj_dir}/data/PPC_Pilot/ATAC-seq
pipe_dir=${prj_dir}/pipeline/ATAC-Seq/sample
script_dir=/projects/PPC/pipeline/ATAC-Seq/work/2021_0928/script

line=`tail -n +$SGE_TASK_ID samples_to_merge.txt | head -1`

l=($line)
id=${l[0]}
r1=${l[1]}
r2=${l[2]}

#rm -rf ${pipe_dir}/${id}

#if [ -d ${pipe_dir}/${id} ]; then exit; fi

#mkdir ${pipe_dir}/${id}

if [ ! -f ${pipe_dir}/${id}/Aligned.out.bam ]; then
cmd="${script_dir}/merge_star.pbs $id $r1 $r2 ${pipe_dir}/${id}"
qsub -N star_${id} $cmd
fi

if [ ! -f ${pipe_dir}/${id}/Aligned.out.filt.bam ]; then 
cmd="${script_dir}/filt.pbs ${pipe_dir}/${id}/Aligned.out.bam ${pipe_dir}/${id}"
qsub -N filt_${id} -hold_jid star_${id} $cmd
fi

if [ ! -f ${pipe_dir}/${id}/plink.genome ]; then 
cmd="${script_dir}/identity.pbs ${pipe_dir}/${id}/Aligned.out.filt.bam ${pipe_dir}/${id}"
qsub -N ident_${id} -hold_jid filt_${id} $cmd
fi

if [ ! -f ${pipe_dir}/${id}/broad_peaks.broadPeak ]; then
cmd="${script_dir}/peak.pbs ${pipe_dir}/${id}/Aligned.out.filt.bam broad ${pipe_dir}/${id}"
qsub -N peak_${id} -hold_jid filt_${id} $cmd
fi

if [ ! -f ${pipe_dir}/${id}/reads_in_promoter.txt ]; then
cmd="${script_dir}/reads_in_promoter.pbs ${pipe_dir}/${id} ${pipe_dir}/${id}"
qsub -N prom_${id} -hold_jid filt_${id} $cmd
fi

if [ ! -f ${pipe_dir}/${id}/broad_peaks.count ]; then 
cmd="${script_dir}/count.pbs ${pipe_dir}/${id}/Aligned.out.filt.bam ${pipe_dir}/${id}/broad_peaks.broadPeak ${pipe_dir}/${id}/broad_peaks.count"
qsub -N count_${id} -hold_jid peak_${id} $cmd
fi

if [ ! -f ${pipe_dir}/${id}/Aligned.out.filt.fs ]; then
cmd="${script_dir}/fsize.pl ${pipe_dir}/${id}/Aligned.out.filt.bam > ${pipe_dir}/${id}/Aligned.out.filt.fs"
echo "date >& 2; echo $cmd >& 2; $cmd; date >& 2" |\
qsub -V -cwd -e logs -o logs -N fsize_${id} -hold_jid filt_${id}
fi

if [ ! -f ${pipe_dir}/${id}/Aligned.out.mdup.q20 ]; then
cmd="samtools view -c -q 20 ${pipe_dir}/${id}/Aligned.out.mdup.bam > ${pipe_dir}/${id}/Aligned.out.mdup.q20"
echo "date >& 2; echo '$cmd' >& 2; $cmd; date >& 2" |\
qsub -V -cwd -e logs -o logs -N q20_${id} -hold_jid filt_${id}
fi

