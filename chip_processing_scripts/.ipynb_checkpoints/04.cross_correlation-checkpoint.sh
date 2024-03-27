#!/bin/bash

#$ -N peaks
#$ -V -cwd
#$ -o logs/cc
#$ -e logs/cc
#$ -pe smp 8

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

set -e

out_dir=$1
upperbound=$2 # 50 for TF, 100 for histone
functions_file=$3

log=${out_dir}/pipeline.log

source $functions_file

if [ ! -d scratch ]; then mkdir scratch; fi 

dir=`mktemp -d -p scratch`

# 1. Make tagAlign for filtered (but not deduped) BAM and subsample for cc analysis (from ENCODE)

# remove low-quality reads
cmd="samtools view -@ 12 -F 1804 -f 2 -q 30 -o ${dir}/Aligned.filt.bam ${out_dir}/Aligned.bam"
print_exec "$cmd" $log

nreads=10000000
ta_file=${dir}/Aligned.merged.tagAlign
subsampled_ta_file=${dir}/Aligned.merged.subsampled.tagAlign
bedtools bamtobed -i ${dir}/Aligned.filt.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > ${ta_file}
grep -v “chrM” ${ta_file} | shuf -n ${nreads} > ${subsampled_ta_file}

# 2. Estimate read length from first 100 reads.
read_len=$(head -n 100 ${ta_file} | awk 'function abs(v) {{return v < 0 ? -v : v}} BEGIN{{sum=0}} {{sum+=abs($3-$2)}} END{{print int(sum/NR)}}')

# 3. Determine exclusion range for fragment length estimation.
# Use a fixed lowerbound at -500.
# Upperbound EXCLUSION_RANGE_MAX is 
#   TF ChIP-seq:  max(read_len + 10, 50)
#   Histone ChIP-seq:  max(read_len + 10, 100)

# lowerbound is fixed at 500 for both
exclusion_range_min=-500
exclusion_range_max=$(python -c "read_len = $read_len; upperbound = $upperbound; print(max(read_len + 10, upperbound))")
#rm -f ${TA_FILE}.tmp

# 4. cross-correlation analysis
nthreads=8
cc_scores_file="${subsampled_ta_file}.cc"
cc_plot_file="${subsampled_ta_file}.cc.plot.pdf"

# cc_score file format

columns="1s/^/Filename\tnumReads\testFragLen\tcorr_estFragLen\tPhantomPeak\tcorr_phantomPeak\targmin_corr\tmin_corr\tphantomPeakCoef\trelPhantomPeakCoef\tQualityTag\n/"

source /frazer01/home/jennifer/.bash_profile
source $functions_file

cmd="Rscript /frazer01/home/jennifer/software/phantompeakqualtools/run_spp.R -c=${subsampled_ta_file} -p=${nthreads} -filtchr=chrM -savp=${cc_plot_file} -out=${cc_scores_file} -x=${exclusion_range_min}:${exclusion_range_max} -rf"
print_exec "$cmd" $log

sed -r 's/,[^\t]+//g' ${cc_scores_file} > ${dir}/temp
mv ${dir}/temp ${cc_scores_file}

sed -i $columns ${cc_scores_file}

cmd="mv ${cc_plot_file} ${cc_scores_file} ${out_dir}/qc"
print_exec "$cmd" $log

cmd="rm -r ${dir}"
print_exec "$cmd" $log
