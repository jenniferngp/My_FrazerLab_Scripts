if [ ! -f ${bed_chip}.gz ] || [ -s ${bed_chip}.gz  ]
then
    
    if [ ! -f ${bam_chip}.nsort.bam ]
    then
        cmd="samtools sort -@ 4 -n -o ${bam_chip}.nsort.bam ${bam_chip}"
        echo $cmd >> $log_file
        echo $cmd >& 2
        eval $cmd
    fi
    
    if [ ! -f ${bam_chip}.nsort.fixed.bam ]
    then
        cmd="samtools fixmate -@ 4 ${bam_chip}.nsort.bam ${bam_chip}.nsort.fixed.bam"
        echo $cmd >> $log_file
        echo $cmd >& 2
        eval $cmd
    fi
    
    if [ ! -f ${bed_chip}.tmp ]
    then
        cmd="bedtools bamtobed -bedpe -mate1 -i ${bam_chip}.nsort.fixed.bam | gzip -nc > ${bed_chip}.tmp"
        echo $cmd >> $log_file
        echo $cmd >& 2
        eval $cmd
    fi

    zcat ${bed_chip}.tmp | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${bed_chip}.gz
fi

if [ -s ${bed_input}.gz ] || [ ! -f ${bed_input}.gz ]
then 

    if [ ! -f ${bam_input}.nsort.bam ]
    then 
        cmd="samtools sort -@ 4 -n -o ${bam_input}.nsort.bam ${bam_input}"
        echo $cmd >> $log_file
        echo $cmd >& 2
        eval $cmd
    fi
    
    if [ ! -f ${bam_input}.nsort.fixed.bam ]
    then
        cmd="samtools fixmate -@ 4 ${bam_input}.nsort.bam ${bam_input}.nsort.fixed.bam"
        echo $cmd >> $log_file
        echo $cmd >& 2
        eval $cmd
    fi

    if [ ! -f ${bed_input}.tmp ]
    then
        cmd="bedtools bamtobed -bedpe -mate1 -i ${bam_input}.nsort.fixed.bam | gzip -nc > ${bed_input}.tmp"
        echo $cmd >> $log_file
        echo $cmd >& 2
        eval $cmd
    fi

    if [ ! -f ${bed_input}.gz ]
    then
        zcat ${bed_input}.tmp | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${bed_input}.gz
    fi
fi

cmd="rm ${bam_chip}.nsort.bam ${bam_input}.nsort.bam ${bed_chip}.tmp ${bed_input}.tmp ${bam_input}.nsort.fixed.bam ${bam_chip}.nsort.fixed.bam"
echo $cmd; eval $cmd