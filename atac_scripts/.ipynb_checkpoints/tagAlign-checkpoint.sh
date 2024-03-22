#!/bin/bash

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

# ATAC-seq pipeline
# Documentation: https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
# Github: https://github.com/ENCODE-DCC/atac-seq-pipeline/tree/master

# ===================
# Create tagAlign file
# ===================
OUT_DIR=$1
PREFIX=Aligned.sorted
FINAL_NMSRT_BAM_PREFIX=${OUT_DIR}/Aligned.sorted.filt.nodup
FINAL_NMSRT_BAM_FILE=${FINAL_NMSRT_BAM_PREFIX}.bam
FINAL_BEDPE_FILE="${FINAL_NMSRT_BAM_PREFIX}.bedpe.gz"
FINAL_TA_FILE="${FINAL_NMSRT_BAM_PREFIX}.tagAlign.gz"
SCRIPT_DIR=/projects/PPC/pipeline/ATAC-Seq/encode/2707259c-f2b0-4ebc-acb1-60aadd0e34d5/scripts

# =================================
# Subsample tagAlign file
# Restrict to one read end per pair for CC analysis
# ================================
NREADS=25000000
SUBSAMPLED_TA_FILE="${PREFIX}.filt.nodup.sample.$((NREADS / 1000000)).MATE1.tagAlign.gz"

zcat ${FINAL_BEDPE_FILE} | grep -v “chrM” | shuf -n ${NREADS} --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f ${FINAL_TA_FILE} | wc -c) -nosalt </dev/zero 2>/dev/null)  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"N","1000",$9}' | gzip -nc > ${SUBSAMPLED_TA_FILE}


# =================================
# Cross-correlation QC
# ================================
CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.cc.qc"
CC_PLOT_FILE="${SUBSAMPLED_TA_FILE}.cc.plot.pdf"
NTHREADS=4

# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag

#Rscript /frazer01/home/jennifer/software/phantompeakqualtools/run_spp.R -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE}

#sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
#mv temp ${CC_SCORES_FILE}

# =================================
# Tn5 shifting
# ================================

SHIFTED_TAG_FILE="${FINAL_NMSRT_BAM_PREFIX}.tn5.tagAlign.gz"
zcat ${FINAL_TA_FILE} | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -nc > $SHIFTED_TAG_FILE

# =================================
# GC bias 
# ================================
export PATH=/software/biobambam2-2.0.95/bin:/software/STAR-2.7.3a/bin/Linux_x86_64:/software/samtools-1.9:/software/sambamba_v0.6.7:/software/rsem-1.2.20:$PATH
REF=/frazer01/home/jennifer/references/Gencode.v34lift37/GRCh37.primary_assembly.genome.fa
GC_BIAS_LOG=${FINAL_NMSRT_BAM_PREFIX}.gc.bias.txt
GC_BIAS_PLOT=${FINAL_NMSRT_BAM_PREFIX}.gc.bias.plot
cmd="java -jar /software/picard-2.20.1/picard.jar CollectGcBiasMetrics R=${REF} I=${FINAL_NMSRT_BAM_FILE} O=${GC_BIAS_LOG} \
USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE \
VERBOSITY=ERROR QUIET=TRUE \
ASSUME_SORTED=FALSE \
CHART=${GC_BIAS_PLOT} S=summary.txt"
echo $cmd; eval $cmd

cmd="python3 ${SCRIPT_DIR}/plot_gc.py $GC_BIAS_LOG"
echo $cmd; eval $cmd

# =================================
# Fragment length statistics
# ================================
INSERT_DATA=${FINAL_NMSRT_BAM_PREFIX}.insertsize.metrics
INSERT_PLOT=${FINAL_NMSRT_BAM_PREFIX}.insertsize.metrics.plot

cmd="java -jar /software/picard-2.20.1/picard.jar CollectInsertSizeMetrics \
INPUT=${FINAL_NMSRT_BAM_FILE} OUTPUT=${INSERT_DATA} H=${INSERT_PLOT} \
VERBOSITY=ERROR QUIET=TRUE \
USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE \
W=1000 STOP_AFTER=5000000"
echo $cmd; eval $cmd

cmd="python3 ${SCRIPT_DIR}/plot_fragment_length.py $INSERT_DATA $FINAL_NMSRT_BAM_PREFIX"
echo $cmd; eval $cmd

cmd="mv $INSERT_DATA $GC_BIAS_LOG *.png *.qc *.flagstat *.tagAlign.gz qc"
echo $cmd; eval $cmd