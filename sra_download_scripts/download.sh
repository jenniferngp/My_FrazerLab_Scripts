

srr=$1
outdir=$2
key=$3

# 1. download the compressed SRA format onto local to avoid download interruption in next step
cmd="prefetch --ngc $key --max-size 420000000000 $srr"
echo $cmd >& 2; eval $cmd

# To look at info
# vdb-dump --info SRR8419492.sra --ngc /frazer01/home/hiroko/dbGaP/prj_20817_D13355.ngc

# More help:
# http://rcs.bu.edu/examples/bioinformatics/sratoolkit/

cd ~/ncbi/public/sra

# 2. extract data from .sra as fastq files
cmd="fasterq-dump --ngc $key --split-files ${srr}.sra"
echo $cmd >& 2; eval $cmd

# 3. rename and move
if [ ! -d ${outdir}/${srr} ]; then mkdir ${outdir}/${srr}; fi
cmd="mv ${srr}_1.fastq ${outdir}/${srr}/${srr}_1.fastq"
echo $cmd >& 2; eval $cmd

cmd="mv ${srr}_2.fastq ${outdir}/${srr}/${srr}_2.fastq"
echo $cmd >& 2; eval $cmd


