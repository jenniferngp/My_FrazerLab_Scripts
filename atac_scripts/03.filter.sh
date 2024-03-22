#!/bin/bash

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

out_dir=$1
log=$2
functions_file=$3

source $functions_file

echo "####### Running 03.filter.sh #######" >> $log

prefix=${out_dir}/Aligned.sorted
raw_bam_file=${prefix}.bam
filt_bam_prefix=${prefix}.filt
filt_bam_file=${filt_bam_prefix}.bam
tmp_filt_bam_prefix=${filt_bam_prefix}.nmsrt
tmp_filt_bam_file=${tmp_filt_bam_prefix}.bam
tmp_filt_fixmate_bam_file=${tmp_filt_bam_prefix}.fixmate.bam
multimapping=4

tmp_filt_bam_file=${filt_bam_prefix}.dupmark.bam
dup_file_qc=${filt_bam_prefix}.dupmark.picard
tmp_filt_bam_file_qc=${filt_bam_prefix}.dupmark.flagstat

final_bam_prefix=${filt_bam_prefix}.nodup
final_bam_file=${final_bam_prefix}.bam # To be stored
final_bam_index_file=${final_bam_file}.bai
final_bam_file_mapstats=${final_bam_prefix}.flagstat # QC file

encode_dir=/frazer01/home/jennifer/software/encode-atac-seq-pipeline/src

# 1. Remove  unmapped, mate unmapped, not primary alignment, reads failing platform
cmd="samtools view -@ 8 -F 524 -f 2 -q 30 -u ${raw_bam_file} | samtools sort -@ 8 -n -o ${tmp_filt_bam_file} -"
print_exec "$cmd" $log

# 2. Take multimappers and randomly assign a location for them
cmd="samtools view -h ${tmp_filt_bam_file} | ${encode_dir}/assign_multimappers.py -k $multimapping --paired-end | samtools fixmate -@ 8 -r - ${tmp_filt_fixmate_bam_file}"
print_exec "$cmd" $log

# 3. Remove reads unmapped, mate unmapped, not primary alignment, failed platform QC, PCR duplicates. Keep reads in proper pair
cmd="samtools view -@ 8 -F 1804 -f 2 -q 30 -u ${tmp_filt_fixmate_bam_file} | samtools sort -@ 8 -o ${filt_bam_file} -"
print_exec "$cmd" $log

cmd="rm ${tmp_filt_fixmate_bam_file} ${tmp_filt_bam_file}"
print_exec "$cmd" $log

# 4. Mark duplicates
cmd="picard MarkDuplicates INPUT=${filt_bam_file} OUTPUT=${tmp_filt_bam_file} METRICS_FILE=${dup_file_qc} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
print_exec "$cmd" $log

# 5. run flagstat to get metrics
cmd="samtools flagstat -@ 10 ${tmp_filt_bam_file} > ${tmp_filt_bam_file_qc}"
print_exec "$cmd" $log

# 6. clean
cmd="mv ${tmp_filt_bam_file} ${filt_bam_file}"
print_exec "$cmd" $log

cmd="mv ${dup_file_qc} ${out_dir}/qc"
print_exec "$cmd" $log

# 7. Remove duplicates
cmd="samtools view -@ 8 -F 1804 -f 2 -b ${filt_bam_file} > ${final_bam_file}"
print_exec "$cmd" $log

# 8. Index Final BAM file
cmd="samtools index -@ 8 ${final_bam_file}"
print_exec "$cmd" $log

# 9. Remove mitochondrial chromosomes
nomito=${filt_bam_prefix}.nodup.nomito.bam
cmd="samtools idxstats ${final_bam_file} | cut -f 1 | grep -v -P "^chrM" | xargs samtools view ${final_bam_file} -@ 8 -b > ${nomito}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="samtools index -@ 8 ${nomito} ${nomito}.bai"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 10. run flagstat to get metrics
cmd="samtools flagstat -@ 8 ${nomito} > ${nomito}.flagstat"
print_exec "$cmd" $log

cmd="mv ${nomito}.flagstat ${out_dir}/qc; rm ${filt_bam_file}"
print_exec "$cmd" $log

