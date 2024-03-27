#$ -N counts
#$ -V -cwd
#$ -pe smp 4
#$ -o logs/counts
#$ -e logs/counts

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

out_dir=$1
saf=$2

log=${out_dir}/logs/counts.log

if [ -f ${out_dir}/Aligned.merged.bam ]
then
	bam=${out_dir}/Aligned.merged.bam
else
	bam=${out_dir}/Aligned.filt.srt.nodup.bam
fi

if [ ! -d ${out_dir}/ref_peaks ]; then mkdir ${out_dir}/ref_peaks; fi

date >& 2
date >> $log

cmd="featureCounts -p -T 4 --donotsort -F SAF -a ${saf} -o ${out_dir}/ref_peaks/ref_peaks.counts $bam"
echo $cmd >> $log
echo $cmd >& 2; eval $cmd
