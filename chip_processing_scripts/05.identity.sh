#$ -N chip-identity
#$ -V
#$ -cwd
#$ -e logs/identity
#$ -o logs/identity
#$ -pe smp 8
#$ -t 6-6:1
#$ -tc 20

set -e

source /frazer01/home/jennifer/.bash_profile
source activate frazer-rna

reference_vcf=/projects/PPC/pipeline/RNA-Seq/work/2023_0804/plink/cardips_common_hg38.vcf.gz 
reference=/reference/public/ucsc/hg38/hg38.fa

out_dir=$1 
functions_file=$2 

source $functions_file

log=${out_dir}/pipeline.log

bam=${out_dir}/Aligned.filt.srt.nodup.bam

dir=`mktemp -d -p scratch` # make scratch directory

echo "$bam target" > ${dir}/sample.txt # make file to rename bamfile with "target" in VCF header

cmd="bcftools mpileup -Ou -f $reference -R $reference_vcf --threads 12 $bam |\
bcftools call --threads 12 -Ou -mv |\
bcftools filter -e 'DP<10'  | bcftools reheader -s ${dir}/sample.txt |\
awk 'BEGIN{OFS=\"\t\";} {if(\$1 !~ /^#/) {\$3=\$1\":\"\$2;} print;}' |\
bgzip -c > ${dir}/call.vcf.gz"
print_exec "$cmd" $log

cmd="tabix -p vcf ${dir}/call.vcf.gz"
print_exec "$cmd" $log

cmd="bcftools merge ${dir}/call.vcf.gz $reference_vcf |\
bcftools view --threads 8 -m2 -M2 -v snps -o ${dir}/merged.vcf.gz -O z"
print_exec "$cmd" $log

cmd="plink --threads 12 --vcf ${dir}/merged.vcf.gz --make-bed --out ${dir}/merged --allow-extra-chr"
print_exec "$cmd" $log

cmd="plink --threads 12 --bfile ${dir}/merged --genome full --out ${dir}/plink --allow-extra-chr"
print_exec "$cmd" $log

cmd="rsync ${dir}/plink.genome ${out_dir}/."
print_exec "$cmd" $log

cmd="rm -rf $dir"
print_exec "$cmd" $log
