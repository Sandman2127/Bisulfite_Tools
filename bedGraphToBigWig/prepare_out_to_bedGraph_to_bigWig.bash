#!/bin/bash

export WD="/mnt/dv/wid/projects1/Zhong-HistoneModification/dsanders/Data/Sequencing_data/05_10_16_JJ_SNPs_Highalt_A_thaliana/DNA_met/lhasa_generated_cytosine_background/for_jbrowse"

#for i in *.out ;

#do

#non stranded bGph separation

#awk '{print $1"\t"$2-1"\t"$2"\t"$5}' $1 | sort -k1,1 -k2,2n > ${1%.out}.bGph

#stranded bGph

#awk '{print $1"\t"$2-1"\t"$2"\t"$3$5}' $i | sort -k1,1 -k2,2n > ${i%.out}.bGph

for i in *.bGph ; 

do

$WD/bedGraphToBigWig $WD/${i%.out}.bGph $WD/tairchrs ${i%.out}.bw ;

done

#done


