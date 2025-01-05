# 2nd run for 2021_0909

prj_dir=/projects/PPC
data_dir=${prj_dir}/data/PPC_Pilot/ATAC-seq
pipe_dir=${prj_dir}/pipeline/ATAC-Seq/sample
script_dir=/projects/PPC/pipeline/ATAC-Seq/work/2021_0928/script

if [ ! -d logs ]; then mkdir logs; fi

# Separate into directories
#igm_dir=/projects/PPC/data/PPC_Pilot/download/igm-storage2.ucsd.edu/210923_A00953_0408_BHVHT7DSX2
#ls ${igm_dir}/*.fastq.gz | while read id; do id=`basename $id`; if [[ ${id} != *Undetermined* ]]; then echo  ${id::-16} >> samples.tmp; fi; done
#sort -u samples.tmp > samples.txt
#rm samples.tmp

#cat samples.txt | while read id; do
#        if [ ! -d ${data_dir}/${id} ]
#	then 
#		mkdir ${data_dir}/${id}
#		mv ${igm_dir}/${id}_* ${data_dir}/${id}
#	else 
#		echo $id exists.
#	fi
#done

# Run pipeline
#nsamples=`wc -l samples.txt`
#qsub -N atac -l short -o logs -e logs -V -cwd -t 1-109:1 ${script_dir}/run_pipeline.sh

# Gather Bam files from Run 1 and 2
#if [ ! -d bam ]; then mkdir bam; fi
#cat all_samples.txt | while read id; do 
#cd bam
#ln -s ${pipe_dir}/${id}/Aligned.out.filt.bam 
#mv Aligned.out.filt.bam ${id}.bam
#cd ..
#done

# Call peaks
#cmd="export PATH=/home/hiroko/Python-3.6.4/bin:$PATH; 
#export LD_LIBRARY_PATH=/home/hiroko/Python-3.6.4/lib:$LD_LIBRARY_PATH;
#macs2 callpeak -f BAMPE -g hs -t bam/*.bam -n merged --outdir `pwd` --nomodel --shift -100 --extsize 200 --broad"
#echo $cmd | qsub -N merged_peaks -o logs -e logs -V -cwd -pe smp 8 -l short

# Call featureCounts
#cmd="${script_dir}/count.pbs bam/*.bam merged.broadPeak `pwd`/merged.broad_peaks.count"
#qsub -N count -hold_jid merged_peaks -l short $cmd

# Gather QC metrics
#id=`cat samples.txt | head -n 1`
#echo -en "id\t" > qc.txt
#script/qc.pl ${pipe_dir}/${id} | head -n 1 >> qc.txt
#echo -en "${id}\t" >> qc.txt
#script/qc.pl ${pipe_dir}/${id} | tail -n +2 >> qc.txt
#cat samples.txt | tail -n +2 | while read id; do
#        echo -en "${id}\t"
#        script/qc.pl ${pipe_dir}/${id} | tail -n +2
#done >> qc.txt

# Change sample names 
#mv S07801_PPC_C5P19_PPC038_ATAC_R01L01S02_S49_L004 S08101_PPC_C3P20_PPC041_ATAC_R01L01S02_S49_L004
#mv S08101_PPC_C3P20_PPC041_ATAC_R01L01S02_S53_L004 S07801_PPC_C5P19_PPC038_ATAC_R01L01S02_S53_L004

#mv S06001_PPC_C1P19_PPC133_ATAC_R02L01S01_S94_L004 P0351_PPC_C2P24_PPC083_ATAC_R02L01S01_S94_L004
#mv P0351_PPC_C2P24_PPC083_ATAC_R01L01S01_S54_L004  S06001_PPC_C1P19_PPC133_ATAC_R02L01S01_S94_L004

#mv P0351_PPC_C2P24_PPC083_ATAC_20190517_S7_L004  S06001_PPC_C1P19_PPC133_ATAC_20190517_S7_L004
#mv S06001_PPC_C1P19_PPC133_ATAC_R02L01S02_S41_L004 P0351_PPC_C2P24_PPC083_ATAC_R02L01S02_S41_L004

# Call narrow peaks
cat samples

cmd="${script_dir}/peak.pbs ${pipe_dir}/${id}/Aligned.out.filt.bam broad ${pipe_dir}/${id}"
qsub -N peak_${id} -hold_jid filt_${id} -l short $cmd
