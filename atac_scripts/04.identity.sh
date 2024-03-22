#$ -V
#$ -cwd
#$ -e logs
#$ -o logs
#$ -pe smp 4

hostname >& 2

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

set -e

out_dir=$1
log=$2
functions_file=$3

plink_dir=${out_dir}/plink
id=`basename $out_dir`
bam=${out_dir}/Aligned.sorted.filt.nodup.nomito.bam
reference_vcf=/projects/PPC/pipeline/RNA-Seq/work/2023_0804/plink/cardips_common_hg38.vcf.gz 
reference=/reference/public/ucsc/hg38/hg38.fa

source $functions_file

echo "####### Running 04.identity.sh #######" >> $log

if [ ! -d ${plink_dir} ]; then mkdir ${plink_dir}; fi
if [ ! -f ${bam}.bai ]; then samtools index ${bam}; fi

# 0. bcftools call will output a vcf with sample header as the bam filename ("Aligned.bam". to make downstream processing easier, we'll rename from "Aligned.bam" to "target"
# to do this, create a text file containing: [old samplename] [new samplename]
echo "${bam} target" > ${plink_dir}/sample.txt

# 1. call genotypes from bam file. outputs a vcf file
cmd="bcftools mpileup --threads 4 -Ou -f $reference -R $reference_vcf $bam |\
bcftools call --threads 4 -Ou -mv |\
bcftools filter -e 'DP<10'  | bcftools reheader -s ${plink_dir}/sample.txt |\
awk 'BEGIN{OFS=\"\t\";} {if(\$1 !~ /^#/) {\$3=\$1\":\"\$2;} print;}' |\
bgzip -c > ${plink_dir}/call.vcf.gz"
print_exec "$cmd" $log

# 2. index vcf file
cmd="tabix -p vcf ${plink_dir}/call.vcf.gz"
print_exec "$cmd" $log

# 3. merge vcf with reference vcf containing wgs calls
cmd="bcftools merge ${plink_dir}/call.vcf.gz $reference_vcf |\
bcftools view --threads 4 -m2 -M2 -v snps -o ${plink_dir}/merged.vcf.gz -O z"
print_exec "$cmd" $log

# 4. convert vcf to plink bed files
cmd="plink --threads 4 --vcf ${plink_dir}/merged.vcf.gz --make-bed --out ${plink_dir}/merged --allow-extra-chr"
print_exec "$cmd" $log

# 5. run plink genome which outputs IBD for each pair of samples in the vcf
cmd="plink --threads 4 --bfile ${plink_dir}/merged --genome full --out ${plink_dir}/plink --allow-extra-chr"
print_exec "$cmd" $log

# 6. clean and remove files
cmd="rm ${plink_dir}/call.vcf.gz ${plink_dir}/call.vcf.gz.tbi ${plink_dir}/sample.txt"
print_exec "$cmd" $log
