#!/bin/bash

#$ -pe smp 8
#$ -V -cwd

# qsub -t 1-[number of uuids]:1 -tc 10 -pe smp 8 -o logs -e logs -V -cwd pipeline.sh [fastq_dir] [uuid_list] [out_dir] [ref_saf] [script_dir]

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

encode_dir=/frazer01/home/jennifer/software/encode-atac-seq-pipeline

fq_dir=$1 # directory containing fastqs. example: /projects/CARDIPS/data/ATAC-Seq/sample
uuid_list=$2 # one-column list of uuids 
uuid=`tail -n +$SGE_TASK_ID $uuid_list | head -1`
out_dir=$3 # output directory. example: /projects/CARDIPS/pipeline/ATAC-Seq/sample_hg38
out_dir=${out_dir}/${uuid}
ref_saf=$4 # saf file containing reference peaks
script_dir=$5 # directory containing atac processing scripts. example: /projects/CARDIPS/analysis/data_organization/scripts/atac_scripts

functions_file=${script_dir}/functions.sh

source $functions_file

log=${out_dir}/pipeline.log

print_exec "mkdir $out_dir $out_dir/detect_adapters $out_dir/logs $out_dir/qc" $log

# 1. Trim adapters 
# Note: ENCODE uses the same adapter seq for both fq1 and fq1; see "detect_most_likely_adapter" function; https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_task_trim_adapter.py
sh ${script_dir}/00.cutadapt.sh $fq_dir $uuid $out_dir $log $functions_file

# 1.1. Index reference genome for alignment (already done! do not run again. takes a very long time.)
#cmd="source activate encode-atac; bowtie2-build ~/references/Gencode.v34lift37/GRCh37.primary_assembly.genome.fa.gz /reference/private/Gencode.v34lift37/bowtie2_index"
#echo $cmd | qsub -N index -V -cwd -o logs/index.out -e logs/index.err -pe smp 8

# 2. Align using BWA-MEM
sh ${script_dir}/01.align.sh $out_dir $log $functions_file

# 4. Calculate % mito reads
sh ${script_dir}/02.mito.sh $out_dir $log $functions_file

# 3. Filter out reads
sh ${script_dir}/03.filter.sh $out_dir $log $functions_file

# 4. Perform sample identity check
log_out=${out_dir}/logs/id.out
og_err=${out_dir}/logs/id.err
qsub -N id -o $log_out -e $log_err -V -cwd ${script_dir}/04.identity.sh $out_dir $log $functions_file

# 5. Generate QC metrics (GC bias and fragment length)
log_out=${out_dir}/logs/qc.out
log_err=${out_dir}/logs/qc.err
qsub -N qc -o $log_out -e $log_err -V -cwd ${script_dir}/05.qc.sh $out_dir $log $functions_file

# 6. Call peaks using MACS2
log_out=${out_dir}/qc/frip.out
log_err=${out_dir}/qc/frip.err

smoothwindow=150 # default on encode
shiftsize=75
gensz=3.0e9 # genome size (hg38)
qsub -hold_jid qc -N peaks -o $log_out -e $log_err -V -cwd ${script_dir}/06.peaks.sh $out_dir $log $functions_file $smoothwindow $shiftsize $gensz

# 7. Calculate TSS Enrichment (TSSE)
cmd="Rscript ${script_dir}/07.tsse_atacqc.R --pipe_dir $out_dir"
print_exec "$cmd" $log

# 8. Reference peak feature counts. Only run this step if you have reference peak set.
if [ $ref_saf != "NA" ]
then
   if [ ! -d ${out_dir}/ref_peaks ]; then mkdir ${out_dir}/ref_peaks; fi
   log_out=${out_dir}/ref_peaks/log_counts.out
   log_err=${out_dir}/ref_peaks/log_counts.err
   qsub -N counts -o $log_out -e $log_err -V -cwd -pe smp 4 ${script_dir}/08.ref_peaks.sh $out_dir $ref_saf $functions_file
else
    echo "Skipping reference peak feature counts"
fi
