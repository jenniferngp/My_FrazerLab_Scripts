source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

out_dir=$1
functions_file=$2
log=${out_dir}/pipeline.log

MAPQ_THRESH=30
SCRIPT_DIR=/projects/CARDIPS/pipeline/ChIP-Seq/encode_script
OFPREFIX=${out_dir}/Aligned
RAW_BAM_FILE=${out_dir}/Aligned.bam
FILT_BAM_PREFIX=${OFPREFIX}.filt.srt
FILT_BAM_FILE=${FILT_BAM_PREFIX}.bam
TMP_FILT_BAM_FILE=${FILT_BAM_PREFIX}.dupmark.bam

source $functions_file

# 1. Remove unmapped, mate unmapped, not in primary alignment, and duplicates (-F 1804). 
# Keep reads in proper pair (-f 2) and have mapping quality > 30 (-q 30).

cmd="samtools view -@ 12 -F 1804 -f 2 -q ${MAPQ_THRESH} -u ${RAW_BAM_FILE} | samtools sort -@ 12 -n - -o ${OFPREFIX}.sorted.bam.tmp" 
print_exec "$cmd" $log

# 2. Fix mate coordinates
cmd="samtools fixmate -@ 12 -r ${OFPREFIX}.sorted.bam.tmp ${OFPREFIX}.fixmate.tmp"
print_exec "$cmd" $log

# 3. Remove unmapped, mate unmapped, not in primary alignment, and duplicates (-F 1804). 
# Keep reads in proper pair (-f 2) and have mapping quality > 30 (-q 30).

cmd="samtools view -@ 12 -F 1804 -f 2 -q ${MAPQ_THRESH} -u ${OFPREFIX}.fixmate.tmp | samtools sort -@ 12 -o ${FILT_BAM_FILE}"
print_exec "$cmd" $log

cmd="rm ${OFPREFIX}.fixmate.tmp ${OFPREFIX}.sorted.bam.tmp"
print_exec "$cmd" $log

# 4. Mark duplicates
DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc"
cmd="picard MarkDuplicates INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
print_exec "$cmd" $log

cmd="mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}"
print_exec "$cmd" $log

# 5. Remove duplicates
FINAL_BAM_PREFIX="${OFPREFIX}.filt.srt.nodup"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" 
FINAL_BAM_INDEX="${FINAL_BAM_PREFIX}.bam.bai" 
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat" 

cmd="samtools view -@ 12 -F 1804 -f 2 -q ${MAPQ_THRESH} -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}"
print_exec "$cmd" $log

# 6. Index
cmd="samtools index -@ 12 ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX}"
print_exec "$cmd" $log

# 7. Flagstat
cmd="samtools flagstat -@ 12 ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}" # recommend using this, faster
print_exec "$cmd" $log

# 8. Clean
cmd="mv ${out_dir}/*.qc ${out_dir}/*.flagstat ${out_dir}/qc"
print_exec "$cmd" $log

cmd="rm $FILT_BAM_FILE"
print_exec "$cmd" $log

# 9. Bigwig
cmd="bamCoverage --numberOfProcessors 12 --bam $FINAL_BAM_FILE --outFileName ${FINAL_BAM_PREFIX}.bw"
print_exec "$cmd" $log
