
source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

id=$1
fq_dir=$2
out_dir=$3

cropped_fastq1=${out_dir}/trimmed/trimmed_R1.fastq
cropped_fastq2=${out_dir}/trimmed/trimmed_R2.fastq
tmp_cropped1=${out_dir}/trimmed/trimmed_R1.tmp
tmp_cropped2=${out_dir}/trimmed/trimmed_R2.tmp

crop_length=151
crop_length_tol=3

CROP=${crop_length}
MINLEN=$((crop_length-crop_length_tol))
NTH=4

# Concatenate fastq files
cmd="cat ${fq_dir}/*R1* > ${out_dir}/trimmed/R1.fastq.gz"
echo $cmd; eval $cmd

cmd="cat ${fq_dir}/*R2* > ${out_dir}/trimmed/R2.fastq.gz"
echo $cmd; eval $cmd

fastq1=${out_dir}/trimmed/R1.fastq.gz
fastq2=${out_dir}/trimmed/R2.fastq.gz

# Trim fastq files
cmd="java -jar /software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads ${NTH} ${fastq1} ${fastq2} ${cropped_fastq1} ${tmp_cropped1} ${cropped_fastq2} ${tmp_cropped2} MINLEN:${MINLEN} CROP:${CROP}"
echo $cmd; eval $cmd

cmd="rm ${tmp_cropped1} ${tmp_cropped2} $fastq1 $fastq2"
echo $cmd; eval $cmd

