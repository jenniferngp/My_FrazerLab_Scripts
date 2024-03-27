#!/bin/bash/

#$ -V -cwd
#$ -pe smp 8

export PATH=/software/cellranger-7.0.1:$PATH
set -e

if [ ! -d scratch ]; then mkdir scratch; fi

# run_cellranger.sh [id] [fastq_dir] [out_dir]
# fastq_dir=/projects/PPC/data/PPC_Pilot
# out_dir=/projects/PPC/pipeline/scRNA-Seq/sample_hg38

dir=`mktemp -d -p scratch`
id=$1 # sample id
data_dir=$2 # directory containing fastq files
out_dir=$3 # output directory
log=${dir}/pipe.log # log file to output error messages

echo $id $data_dir $out_dir

# 1. Transfer fastqs to scratch directory
date >& 2
cmd="mkdir ${dir}/fastqs; rsync -r ${data_dir}/${id} ${dir}/fastqs"
echo $cmd >> $log
echo $cmd >& 2; eval $cmd

# 2. Get reference
# to download the reference: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation?src=social&lss=facebook&cnm=soc-fb-ra_g-program&cid=NULL
date >& 2
ref_dir=/frazer01/home/jennifer/references/CellRanger/refdata-gex-GRCh38-2020-A/refdata-gex-GRCh38-2020-A

# 3. Run Cellranger
date >& 2
cmd="cellranger count --id=${id} --transcriptome=${ref_dir} --fastqs=${dir}/fastqs --sample=${id} --localcores=8 --localmem=64"
echo $cmd >> $log
echo $cmd >& 2; eval $cmd

# 4. Move output to final directory
date >& 2
cmd="mv ${id} ${out_dir}/"
echo $cmd >> $log
echo $cmd  >& 2; eval $cmd

# 5. Remove scratch
date >& 2
cmd="mv $log ${out_dir}/${id} rm -r ${dir}"
echo $cmd >> $log
echo $cmd >& 2; eval $cmd

date >& 2
echo Done!
