#$ -N rna
#$ -pe smp 16
#$ -V -cwd

set -e

hostname >& 2

source /frazer01/home/jennifer/.bash_profile
source activate frazer-rna
source /projects/CARDIPS/analysis/data_organization/scripts/functions.sh

export PATH=/software/biobambam2-2.0.95/bin:/software/sambamba_v0.6.7:$PATH
reference_star=/reference/private/STAR/hg38_gencode44/reference
reference_rsem=/reference/private/RSEM/hg38_gencode44/reference/hg38_gencode44

if [ $# -ne 3 ]
then
    echo "Exiting. Need input: sh merge_pipeline.sh <fastq_dir> <uuid_list> <out_dir>" >& 2
    exit 1
fi

fq_dir=$1 # fastq directory
uuid=`tail -n +$SGE_TASK_ID $2 | head -1` 
out_dir=$3 # output directory
out_dir=${out_dir}/${uuid}

log=${out_dir}/pipeline.log

if [ ! -d $out_dir ]; then mkdir $out_dir; fi

echo Processing .. ${uuid} >& 2
echo Fastq directory.. ${fq_dir} >& 2
echo Output directory.. ${out_dir} >& 2
echo Working directory.. ${dir} >& 2

# 0. Make scratch directory
if [ ! -d scratch ]; then mkdir scratch; fi
dir=`mktemp -d -p scratch`
if [ ! -d ${dir} ]; then exit 1; fi

# 0. Copy fastqs to scratch directory
cmd="rsync -r ${fq_dir}/${uuid} $dir"
print_exec "$cmd" $log

# 1. Concatenate R1 fastqs from different sequencing runs
cmd="zcat ${dir}/${uuid}/*R1*.fastq.gz > ${dir}/R1.fastq"
print_exec "$cmd" $log

# 2. Concatenate R2 fastqs from different sequencing runs
cmd="zcat ${dir}/${uuid}/*R2*.fastq.gz > ${dir}/R2.fastq"
print_exec "$cmd" $log

# 3. Clean
cmd="rm ${dir}/${uuid}/*.fastq.gz"
print_exec "$cmd" $log

# 4. Align using STAR
# good STAR manual: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

in_files=( ${dir}/R1.fastq ${dir}/R2.fastq )
echo fastq input.. ${in_files[@]} >& 2

cmd="STAR \
--runThreadN 16 \
--genomeDir ${reference_star} \
--genomeLoad NoSharedMemory \
--readFilesIn ${in_files[@]} \
--outSAMattributes All \
--outSAMunmapped Within \
--outSAMattrRGline ID:1 PL:ILLUMINA PU:CARDIPS LB:${uuid} SM:${uuid} \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ${dir}/ \
--quantMode TranscriptomeSAM"
print_exec "$cmd" $log

# 5. Remove fastqs to clear space
cmd="rm ${dir}/R1.fastq ${dir}/R2.fastq"
echo $cmd >& 2; eval $cmd

# 4. Sort bam 
cmd="samtools sort -m 2G -n -o ${dir}/Aligned.out.namesorted.bam -@ 16 ${dir}/Aligned.out.bam"
print_exec "$cmd" $log

# 5. Fill in mate coordinates
cmd="samtools fixmate -@ 16 -m ${dir}/Aligned.out.namesorted.bam ${dir}/Aligned.out.namesorted.fixmate.bam"
print_exec "$cmd" $log

# 6. Sort bam
cmd="samtools sort -m 2G -o ${dir}/Aligned.out.sorted.bam -@ 16 ${dir}/Aligned.out.namesorted.fixmate.bam"
print_exec "$cmd" $log

# 7. Remove bams to clear space
cmd="rm ${dir}/Aligned.out.namesorted.bam ${dir}/Aligned.out.namesorted.fixmate.bam"
print_exec "$cmd" $log

# 8. Index bam
cmd="samtools index -@ 16 ${dir}/Aligned.out.sorted.bam ${dir}/Aligned.out.sorted.bam.bai"
print_exec "$cmd" $log

# 9. Mark duplicates
# Why we don't remove duplicates in RNA: https://dnatech.genomecenter.ucdavis.edu/faqs/should-i-remove-pcr-duplicates-from-my-rna-seq-data/
cmd="samtools markdup -@ 16 -s -T ${dir}/temp ${dir}/Aligned.out.sorted.bam ${dir}/Aligned.out.sorted.mdup.bam"
print_exec "$cmd" $log

# 10. Clean 
cmd="rm ${dir}/Aligned.out.sorted.bam ${dir}/Aligned.out.sorted.bam.bai"
print_exec "$cmd" $log

# 8. Index bam
cmd="samtools index -@ 16 ${dir}/Aligned.out.sorted.mdup.bam ${dir}/Aligned.out.sorted.mdup.bam.bai"
print_exec "$cmd" $log

# 9. Run flagstat
cmd="samtools flagstat -@ 16 ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.bam.flagstat"
print_exec "$cmd" $log

# 10. Run idxstats
cmd="samtools idxstats -@ 16 ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.bam.idxstats"
print_exec "$cmd" $log

# 10. Run samtools stats
cmd="samtools stats -@ 16 ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.bam.stats"
print_exec "$cmd" $log

# 11. Run RSEM (long runtime)
cmd="rsem-calculate-expression \
--bam \
--num-threads 16 \
--no-bam-output \
--seed 3272015 \
--estimate-rspd \
--forward-prob 0 \
--paired-end"
cmd=$cmd" ${dir}/Aligned.toTranscriptome.out.bam ${reference_rsem} ${dir}/rsem"
print_exec "$cmd" $log

# 12. Reformat header for picard
cmd="samtools view -H ${dir}/Aligned.out.sorted.mdup.bam | sed 's,^@RG.*,@RG\tID:None\tSM:None\tLB:None\tPL:Illumina,g' | samtools reheader - ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.reheader.bam"
print_exec "$cmd" $log

# 13. Run Picard's CollectRnaSeqMetrics
refFlat=/reference/public/ucsc/hg38/refFlat.txt.gz

cmd="picard CollectRnaSeqMetrics \
I=${dir}/Aligned.out.sorted.mdup.reheader.bam \
O=${dir}/picard.RNA_Metrics \
VALIDATION_STRINGENCY=SILENT \
REF_FLAT=$refFlat \
STRAND=NONE"
print_exec "$cmd" $log

# 14. Clean
cmd="rm -r ${dir}/Aligned.out.bam ${dir}/Aligned.out.sorted.mdup.reheader.bam ${dir}/Aligned.toTranscriptome.out.bam"
print_exec "$cmd" $log

cmd="mv ${dir}/* ${out_dir}/"
print_exec "$cmd" $log

cmd="rm -r ${dir}"
print_exec "$cmd" $log

