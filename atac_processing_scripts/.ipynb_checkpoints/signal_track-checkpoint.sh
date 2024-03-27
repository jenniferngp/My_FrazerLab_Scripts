#!/bin/bash

out_dir=$1

prefix=${out_dir}/peaks/narrow
fc_bedgraph="$prefix.fc.signal.bedgraph"
fc_bedgraph_srt="$prefix.fc.signal.srt.bedgraph"    
fc_bigwig="$prefix.fc.signal.bigwig"
pval_bedgraph="$prefix.pval.signal.bedgraph"
pval_bedgraph_srt="$prefix.pval.signal.srt.bedgraph"    
pval_bigwig="$prefix.pval.signal.bigwig"

# Make signal track
cmd="macs2 bdgcmp -t $prefix_treat_pileup.bdg -c $prefix_control_lambda.bdg --o-prefix $prefix -m FE"
echo $cmd; eval $cmd

# clip to confine in human genome
chrsz=/reference/public/hg19/hg19.size.txt
cmd="slopBed -i ${prefix}_FE.bdg -g $chrsz -b 0 | bedClip stdin $chrsz ${fc_bedgraph}; rm $prefix_FE.bdg"
echo $cmd; eval $cmd

# Make bedWig
cmd="sort -k1,1 -k2,2n $fc_bedgraph > $fc_bedgraph_srt; bedGraphToBigWig $fc_bedgraph_srt $chrsz $fc_bigwig"
echo $cmd; eval $cmd

# sval counts the number of tags per million in the (compressed) BED file
sval=$(wc -l <(zcat -f "$tag") | awk '{printf "%f", $1/1000000}')

cmd="macs2 bdgcmp -t ${prefix}_treat_pileup.bdg -c ${prefix}_control_lambda.bdg --o-prefix ${prefix} -m ppois -S ${sval}"
echo $cmd; eval $cmd

cmd="slopBed -i ${prefix}_ppois.bdg -g $chrsz -b 0 | bedClip stdin $chrsz $pval_bedgraph"
echo $cmd; eval $cmd

cmd="sort -k1,1 -k2,2n $pval_bedgraph > $pval_bedgraph_srt; bedGraphToBigWig $pval_bedgraph_srt $chrsz $pval_bigwig"
echo $cmd; eval $cmd

#rm -f $pval_bedgraph $pval_bedgraph_srt
#rm -f "$prefix"_treat_pileup.bdg "$prefix"_control_lambda.bdg
#rm -f $fc_bedgraph $fc_bedgraph_srt $prefix_FE.bdg
#rm -f "$prefix"_ppois.bdg


