#!/bin/bash

#Environment Variables:
WD=`pwd -P`
BS_PROGS="/home/BIOTECH/dmsanders/progs/Bisulfite_Tools"
PLOTTER="$BS_PROGS/Full_genome_methylation_binner"
BSMAP="$BS_PROGS/bsmap-2.90/bsmap"
BGPH2BW="$BS_PROGS/bedGraphToBigWig/bedGraphToBigWig"
HCDMR_MAIN="$BS_PROGS/hcDMR_caller/Main"
HCDMR_DATA="$BS_PROGS/hcDMR_caller/Ref_data/"
FASTP="$HOME/progs/fastp/fastp"
CPUCORES=5                          # How many cores do you want to utilize per run?
CG_MET_DIFF=0.4                     # 40% methylation difference
CHG_MET_DIFF=0.2                    # 20% methylation difference
CHH_MET_DIFF=0.1                    # 10% methylation difference
TRIM="TRUE"                         # engage Trimming via fastp or not

#Arguments:
PRIMARY_INPUT=$1
INPUT="$(basename $PRIMARY_INPUT)"  # This allows me to call from a directory structure keeping the data name separate
REF="$2"                            # Reference genome.fa
METHOD="$3"                         # ALL MAP_ONLY CALC_MET TRACKS_ONLY HCDMR_ONLY BIN_ONLY
CONTEXT="$4"                        # ALL_CONTEXTS CG CHG CHH the default is all if no $4, pass CG as $4 with HCDMR only to only call CG DMRs
if [ "$CONTEXT" != "" ]; then
    CONTEXT=$4
    echo "If method is ALL or HCDMR_ONLY calling DMRs on only $CONTEXT"
else
    CONTEXT="ALL_CONTEXTS"
    echo "If method is ALL or HCDMR_ONLY calling DMRs on $CONTEXT"
fi

export PATH HOME WD BS_PROGS PLOTTER BSMAP BGPH2BW CONTEXT

### Build initial directory structure

echo "User decided on the data analysis method: $METHOD"

if [ -e "$WD/${INPUT%.fq.gz}_analysis" ] ; then echo "Initial directory previously made, doing nothing" ; else echo "Making primary analysis directory $WD/${INPUT%.fq.gz}_analysis" ; mkdir ./"${INPUT%.fq.gz}_analysis" ; fi 


### Fastqc all each file
#fastqc $PRIMARY_INPUT --outdir $WD/fastqcs

if ([ "$TRIM" == "TRUE" ] && [ $METHOD == "ALL" ]) || ([ "$TRIM" == "TRUE" ] && [ $METHOD == "MAP_ONLY" ]); then
    TRIMMED_FASTQ="$WD/${INPUT%.fq.gz}_analysis/trimmed_fastq/${INPUT%.fq.gz}.trimmed.fq.gz"
    
    if [ -e "$WD/${INPUT%.fq.gz}_analysis/trimmed_fastq" ]; 
    then 
        echo "$WD/${INPUT%.fq.gz}_analysis/trimmed_fastq already exists, proceeding to trim";
        $FASTP --thread $CPUCORES --in1 $PRIMARY_INPUT --out1 $TRIMMED_FASTQ
        mv ./fastp.json ./fastp.html "$WD/${INPUT%.fq.gz}_analysis/trimmed_fastq"  
    else
        echo "$WD/${INPUT%.fq.gz}_analysis/trimmed_fastq does not exist, making directory and proceeding to trim";
        mkdir -p "$WD/${INPUT%.fq.gz}_analysis/trimmed_fastq" ;
        $FASTP --thread $CPUCORES --in1 $PRIMARY_INPUT --out1 $TRIMMED_FASTQ
        mv ./fastp.json ./fastp.html "$WD/${INPUT%.fq.gz}_analysis/trimmed_fastq"   
    fi
else
    echo "Trimming via FASTP not requested, proceeding to alignment"
fi

if [ "$METHOD" == "ALL" ] || [ "$METHOD" == "MAP_ONLY" ]; 
then
    echo "Bisulfite Mapping engaged with $CPUCORES"
    if [ "$TRIM" == "TRUE" ];
    then
        $BSMAP -f 5 -q 33 -n 0 -a $TRIMMED_FASTQ -d $REF -o ${INPUT%.fq.gz}.bam -p $CPUCORES -v 0.08 ;
    else
        $BSMAP -f 5 -q 33 -n 0 -a $PRIMARY_INPUT -d $REF -o ${INPUT%.fq.gz}.bam -p $CPUCORES -v 0.08 ;
    fi

        # v = 0.12 or 12 % of the 100 bp read can be mismatched, 0.08 = default 8% incorrect
        # -f filter 5 n or less per read
        # -L map first 70bp
        # -q 33 trim anything not 33 quality or better
        # -v 0.08 
        # -n  [0,1]   set mapping strand information. default: 0
        #                   -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+), 
        #                   for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --.
        #                   -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, -- 

    # Base is python 2.7.15
    source activate base 

    python $BS_PROGS/bsmap-2.90/methratio.py -z -u -r -s $BS_PROGS/bsmap-2.90/samtools/ -o ${INPUT%.fq.gz}.outy -d $REF ${INPUT%.fq.gz}.bam ; 
        # Options:
        #         -h, --help            show this help message and exit
        #         -o FILE, --out=FILE   output methylation ratio file name. [default: STDOUT]
        #         -O FILE, --alignment-copy=FILE
        #                                 save a copy of input alignment for BSMAP pipe input.
        #                                 (in BAM format) [default: none]
        #         -w FILE, --wig=FILE   output methylation ratio wiggle file. [default: none]
        #         -b BIN, --wig-bin=BIN
        #                                 wiggle file bin size. [default: 25]
        #         -d FILE, --ref=FILE   reference genome fasta file. (required)
        #         -c CHR, --chr=CHR     process only specified chromosomes, separated by ','.
        #                                 [default: all] example: --chroms=chr1,chr2
        #         -s PATH, --sam-path=PATH
        #                                 path to samtools. [default: none]
        #         -u, --unique          process only unique mappings/pairs.
        #         -p, --pair            process only properly paired mappings.
        #         -z, --zero-meth       report loci with zero methylation ratios.
        #                                 (depreciated, -z will be always enabled)
        #         -q, --quiet           don't print progress on stderr.
        #         -r, --remove-duplicate
        #                                 remove duplicated reads.

    mkdir -p "$WD/${INPUT%.fq.gz}_analysis/Bisulfite_alignment" ;
    mv ${INPUT%.fq.gz}.bam "$WD/${INPUT%.fq.gz}_analysis/Bisulfite_alignment" 
 
else 
    echo "Mapping not engaged proceeding to Genome Binning, to engage mapping use method ALL or MAP_ONLY" ;
fi

if [ "$METHOD" == "ALL" ] || [ "$METHOD" == "CALC_MET" ]; 
then 
    # Calc total methylation levels
    echo "Calculating Genome Wide Methylation Levels in All Contexts (CG, CHG, CHH)"

    _CG=$(awk '$4=="CG"{MC+=$7;TO+=$8;}END{print MC/TO}' ${INPUT%.fq.gz}.outy) ;
    echo "CG mC level (% of total CG):" >> ${INPUT%.fq.gz}_Genome_Met.txt ;
    echo $_CG*100.0 | bc -l >> ${INPUT%.fq.gz}_Genome_Met.txt ;

    _CHG=$(awk '$4=="CHG"{MC+=$7;TO+=$8;}END{print MC/TO}' ${INPUT%.fq.gz}.outy) ;
    echo "CHG mC level (% of total CHG):" >> ${INPUT%.fq.gz}_Genome_Met.txt ;
    echo $_CHG*100.0 | bc -l >> ${INPUT%.fq.gz}_Genome_Met.txt ; 

    _CHH=$(awk '$4=="CHH"{MC+=$7;TO+=$8;}END{print MC/TO}' ${INPUT%.fq.gz}.outy) ;
    echo "CHH mC level (% of total CHH):" >> ${INPUT%.fq.gz}_Genome_Met.txt ;
    echo $_CHH*100.0 | bc -l >> ${INPUT%.fq.gz}_Genome_Met.txt ;

    # Calculate Chloroplast methylation levels to find methylation conversion efficiency
    echo "Calculating Chloroplast Methylation Levels in All Contexts (CG, CHG, CHH, ALL) to estimate bisulfite conversion efficiecy"

    C_CG=$(awk '$1=="chrC" && $4=="CG"{MC+=$7;TO+=$8;}END{print MC/TO}' ${INPUT%.fq.gz}.outy) ;
    echo "CG mC level (% of total CG):" >> ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo $C_CG*100.0 | bc -l >>  ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo "Unmethylated CG Conversion Efficiency equals:" >>  ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo "100.0-($C_CG*100)" | bc -l >>  ${INPUT%.fq.gz}_chrC_Met.txt ;

    C_CHG=$(awk '$1=="chrC" && $4=="CHG"{MC+=$7;TO+=$8;}END{print MC/TO}' ${INPUT%.fq.gz}.outy) ;
    echo "CHG mC level (% of total CHG):" >> ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo $C_CHG*100.0 | bc -l >>  ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo "Unmethylated CHG Conversion Efficiency equals:" >>  ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo "100.0-($C_CHG*100)" | bc -l >>  ${INPUT%.fq.gz}_chrC_Met.txt ;

    C_CHH=$(awk '$1=="chrC" && $4=="CHH"{MC+=$7;TO+=$8;}END{print MC/TO}' ${INPUT%.fq.gz}.outy) ;
    echo "CHH mC level (% of total CHH):" >> ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo $C_CHH*100.0 | bc -l >>  ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo "Unmethylated CHH Conversion Efficiency equals:" >>  ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo "100.0-($C_CHH*100)" | bc -l >>  ${INPUT%.fq.gz}_chrC_Met.txt ;

    C_ALL=$(awk '$1=="chrC"{MC+=$7;TO+=$8;}END{print MC/TO}' ${INPUT%.fq.gz}.outy) ;
    echo "All mC level (% of total C in all contexts):" >> ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo $C_ALL*100.0 | bc -l >>  ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo "Unmethylated Total C Conversion Efficiency equals:" >>  ${INPUT%.fq.gz}_chrC_Met.txt ;
    echo "100.0-($C_ALL*100)" | bc -l >>  ${INPUT%.fq.gz}_chrC_Met.txt ;

else
    echo "Genome Wide methylation levels not calculated, to engage use method ALL or CALC_MET"
fi



if [ "$METHOD" == "ALL" ] || [ "$METHOD" == "TRACKS_ONLY" ] ; 
then
    echo "Building Tracks" ;

    # Non stranded bGph separation
    awk '$4=="CHH"{print $1"\t"$2-1"\t"$2"\t"$5}' ${INPUT%.fq.gz}.outy | awk '$2!=-1{print $0}' | sort -k1,1 -k2,2n > ${INPUT%.fq.gz}.CHH.NS.bGph
    awk '$4=="CHG"{print $1"\t"$2-1"\t"$2"\t"$5}' ${INPUT%.fq.gz}.outy | awk '$2!=-1{print $0}' | sort -k1,1 -k2,2n > ${INPUT%.fq.gz}.CHG.NS.bGph
    awk '$4=="CG"{print $1"\t"$2-1"\t"$2"\t"$5}' ${INPUT%.fq.gz}.outy | awk '$2!=-1{print $0}' | sort -k1,1 -k2,2n > ${INPUT%.fq.gz}.CG.NS.bGph

    $BGPH2BW ${INPUT%.fq.gz}.CHH.NS.bGph $BS_PROGS/bedGraphToBigWig/tairchrs ${INPUT%.fq.gz}.CHH.NS.bw ;
    $BGPH2BW ${INPUT%.fq.gz}.CHG.NS.bGph $BS_PROGS/bedGraphToBigWig/tairchrs ${INPUT%.fq.gz}.CHG.NS.bw ;
    $BGPH2BW ${INPUT%.fq.gz}.CG.NS.bGph $BS_PROGS/bedGraphToBigWig/tairchrs ${INPUT%.fq.gz}.CG.NS.bw ;

    # Stranded bGph
    awk '$4=="CHH"{print $1"\t"$2-1"\t"$2"\t"$3$5}' ${INPUT%.fq.gz}.outy | awk '$2!=-1{print $0}' | sort -k1,1 -k2,2n > ${INPUT%.fq.gz}.CHH.bGph ;
    awk '$4=="CHG"{print $1"\t"$2-1"\t"$2"\t"$3$5}' ${INPUT%.fq.gz}.outy | awk '$2!=-1{print $0}' | sort -k1,1 -k2,2n > ${INPUT%.fq.gz}.CHG.bGph ;
    awk '$4=="CG"{print $1"\t"$2-1"\t"$2"\t"$3$5}' ${INPUT%.fq.gz}.outy | awk '$2!=-1{print $0}' | sort -k1,1 -k2,2n > ${INPUT%.fq.gz}.CG.bGph ;

    $BS_PROGS/bedGraphToBigWig/bedGraphToBigWig ${INPUT%.fq.gz}.CHH.bGph $BS_PROGS/bedGraphToBigWig/tairchrs ${INPUT%.fq.gz}.CHH.bw ;
    $BS_PROGS/bedGraphToBigWig/bedGraphToBigWig ${INPUT%.fq.gz}.CHG.bGph $BS_PROGS/bedGraphToBigWig/tairchrs ${INPUT%.fq.gz}.CHG.bw ;
    $BS_PROGS/bedGraphToBigWig/bedGraphToBigWig ${INPUT%.fq.gz}.CG.bGph $BS_PROGS/bedGraphToBigWig/tairchrs ${INPUT%.fq.gz}.CG.bw ;

        #tairchrs  (i.e your original genome must look like this in the fasta >chr1  ... >chr2)
        # chr1	30427671
        # chr2	19698289
        # chr3	23459830
        # chr4	18585056
        # chr5	26975502
        # chrC	154478
        # chrM	366924

    if [ -s ${INPUT%.fq.gz}.CHH.bw ] && [ -s ${INPUT%.fq.gz}.CHG.bw ] && [ -s ${INPUT%.fq.gz}.CHH.bw ] && [ -s ${INPUT%.fq.gz}.CHH.NS.bw ] && [ -s ${INPUT%.fq.gz}.CHG.NS.bw ] && [ -s ${INPUT%.fq.gz}.CHH.NS.bw ] ; 
    then
        rm ${INPUT%.fq.gz}.CG.bGph ${INPUT%.fq.gz}.CHG.bGph ${INPUT%.fq.gz}.CHH.bGph ${INPUT%.fq.gz}.CG.NS.bGph ${INPUT%.fq.gz}.CHG.NS.bGph ${INPUT%.fq.gz}.CHH.NS.bGph ;
    fi
	
    mkdir -p $WD/"${INPUT%.fq.gz}_analysis"/bigWigs ;
    mv ${INPUT%.fq.gz}*.bw $WD/"${INPUT%.fq.gz}_analysis"/bigWigs 
else
    echo "Track Building NOT engaged, to engage use ALL or TRACKS_ONLY " ;
fi

if [ "$METHOD" == "ALL" ] || [ "$METHOD" == "HCDMR_ONLY" ];
then

    ## Run hcDMR caller
    echo "Running methratio_alt.py from HCDMR" ;

    source activate base  # gets proper python2.7 environment 

    ### Step 1 - generate methratio file from aligned bam file:

    python $HCDMR_MAIN/methratio_alt.py --Steve --sam-path=$BS_PROGS/bsmap-2.90/samtools --ref=$REF --out=${INPUT%.fq.gz}.out -u -z -r "$WD/${INPUT%.fq.gz}_analysis/Bisulfite_alignment/${INPUT%.fq.gz}.bam" ;
        # Options:
        #         -h, --help            show this help message and exit
        #         -o FILE, --out=FILE   output file name. (required)
        #         -d FILE, --ref=FILE   reference genome fasta file. (required)
        #         -c CHR, --chr=CHR     process only specified chromosomes, separated by ','.
        #                                 [default: all] example: --chroms=chr1,chr2
        #         -s PATH, --sam-path=PATH
        #                                 path to samtools. [default: none]
        #         -u, --unique          process only unique mappings/pairs.
        #         -p, --pair            process only properly paired mappings.
        #         -z, --zero-meth       report loci with zero methylation ratios.
        #         -q, --quiet           don't print progress on stderr.
        #         -r, --remove-duplicate
        #                                 remove duplicated reads.
        #         -t N, --trim-fillin=N
        #                                 trim N end-repairing fill-in nucleotides. [default: 2]
        #         -g, --combine-CpG     combine CpG methylaion ratios on both strands.
        #         -m FOLD, --min-depth=FOLD
        #                                 report loci with sequencing depth>=FOLD. [default: 1]
        #         -n, --no-header       don't print a header line
        #         -S, --Steve           activate the Steve filter!!!

    ### Step 2 - count the C and CT count at every position in the genome

    echo "Counting C's and CT's @ every genome position via HCDMR's BSmap_to_cytosine.pl" ;

    perl $HCDMR_MAIN/BSmap_to_cytosine.pl --input_file ${INPUT%.fq.gz}.out.gz --reference_cytosine $HCDMR_DATA/TAIR10_v2.cytosine.gz ;

    ### Step 3 - Bin the genome into 100bp bins.

    echo "Binning genome into 100bp bins via HCDMR's Cytosine_to_100bp.pl" ;

    perl $HCDMR_MAIN/Cytosine_to_100bp.pl ${INPUT%.fq.gz}.out.cytosine.gz ;

    ### Step 4 - Call DMRS: 

    ### Call CHH DMRS :
    if [ "$CONTEXT" == "ALL_CONTEXTS" ] || [ "$CONTEXT" == "CHH" ] ;
    then
        echo "Calling CHH methylation differences with HCDMR" ;

        perl $HCDMR_MAIN/hcDMR_caller_DSmod2.pl -ref $HCDMR_DATA/CHH.100.54WT.Ref.txt.gz -input ${INPUT%.fq.gz}.out.CHH.100.gz -dif $CHH_MET_DIFF -n 33 ;
Col_pUV_BS_seq.out.DMR
        mv ${INPUT%.fq.gz}.out.DMR ${INPUT%.fq.gz}.CHH.DMR ;

        awk '$7>=0.10 && $2!="begin" && $5<$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${INPUT%.fq.gz}.CHH.DMR | grep "hypo" > ${INPUT%.fq.gz}.CHH.hypo.DMR ;
        awk '$7>=0.10 && $2!="begin" && $5>$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${INPUT%.fq.gz}.CHH.DMR | grep "hyper"  > ${INPUT%.fq.gz}.CHH.hyper.DMR ;
    else
        echo "NOT calling CHH methylation differences with HCDMR, CONTEXT is $CONTEXT" ;
    fi 
     ### Call CHG DMRS :
    if [ "$CONTEXT" == "ALL_CONTEXTS" ] || [ "$CONTEXT" == "CHG" ] ;
    then
        echo "Calling CHG methylation differences with HCDMR" ;
        
        perl $HCDMR_MAIN/hcDMR_caller_DSmod2.pl -ref $HCDMR_DATA/CHG.100.54WT.Ref.txt.gz -input ${INPUT%.fq.gz}.out.CHG.100.gz -dif $CHG_MET_DIFF -n 33 ;
        mv ${INPUT%.fq.gz}.out.DMR ${INPUT%.fq.gz}.CHG.DMR  ;
        
        awk '$7>=0.2 && $2!="begin" && $5<$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${INPUT%.fq.gz}.CHG.DMR | grep "hypo" > ${INPUT%.fq.gz}.CHG.hypo.DMR  ;
        awk '$7>=0.2 && $2!="begin" && $5>$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${INPUT%.fq.gz}.CHG.DMR | grep "hyper" > ${INPUT%.fq.gz}.CHG.hyper.DMR  ;
    else
        echo "NOT calling CHG methylation differences with HCDMR, CONTEXT is $CONTEXT" ;
    fi

    ### Call CG DMRS :
    if [ "$CONTEXT" == "ALL_CONTEXTS" ] || [ "$CONTEXT" == "CG" ] ;
    then
        echo "Calling CG methylation differences with HCDMR" ; 
        
        perl $HCDMR_MAIN/hcDMR_caller_DSmod2.pl -ref $HCDMR_DATA/CG.100.54WT.Ref.txt.gz -input ${INPUT%.fq.gz}.out.CG.100.gz -dif $CG_MET_DIFF -n 33 ;
        mv ${INPUT%.fq.gz}.out.DMR ${INPUT%.fq.gz}.CG.DMR ;
        
        awk '$7>=0.4 && $2!="begin" && $5<$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${INPUT%.fq.gz}.CG.DMR | grep "hypo" > ${INPUT%.fq.gz}.CG.hypo.DMR  ;
        awk '$7>=0.4 && $2!="begin" && $5>$6{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t""+"}' ${INPUT%.fq.gz}.CG.DMR | grep "hyper" > ${INPUT%.fq.gz}.CG.hyper.DMR ;
    
    else
        echo "NOT calling CHG methylation differences with HCDMR, CONTEXT is $CONTEXT" ;
    fi 

    ### Summarize DMRs
    echo "Summarizing all hypo and hyper called DMRs in ${INPUT%.fq.gz}_DMR_Count_Summary.DMR" ;
    for _FILE in *hypo.DMR *hyper.DMR ; 
    do 
        echo "${_FILE%.DMR} DMR count:" >> ${INPUT%.fq.gz}_DMR_Count_Summary.DMR ;
        FILECOUNT=$(wc -l $_FILE | awk '{print $1}')
        echo "$FILECOUNT" >> ${INPUT%.fq.gz}_DMR_Count_Summary.DMR ;
        FILECOUNT=""
    done    

    ### Making bw tracks of DMR files:

    for i in ${INPUT%.fq.gz}.CHH.hyper.DMR ${INPUT%.fq.gz}.CHH.hypo.DMR ${INPUT%.fq.gz}.CHG.hyper.DMR ${INPUT%.fq.gz}.CHG.hypo.DMR ${INPUT%.fq.gz}.CG.hyper.DMR ${INPUT%.fq.gz}.CG.hypo.DMR ; 
    do
        awk '{print $1"\t"$2"\t"$3"\t"$6$5}' $i | sort -k1,1 -k2,2n > ${i%.DMR}.bGph ;
        $BS_PROGS/bedGraphToBigWig/bedGraphToBigWig ${i%.DMR}.bGph $BS_PROGS/bedGraphToBigWig/tairchrs ${i%.DMR}.bw ;
        rm ${i%.DMR}.bGph ;
    done
else
    echo "hcDMR calling not engaged, you chose the method $METHOD " ;
fi


if [ "$METHOD" == "ALL" ];
then
    echo "Final cleanup: moving DMR files, transformed DMR bigWigs, outy files and gzipped intermediates to $WD/${INPUT%.fq.gz}_analysis/hcDMRs"
    mkdir -p $WD/"${INPUT%.fq.gz}_analysis"/hcDMRs ;
    mv ${INPUT%.fq.gz}*hyp*.bw ${INPUT%.fq.gz}*.DMR ${INPUT%.fq.gz}*.gz "$WD/${INPUT%.fq.gz}_analysis/hcDMRs"
    mv ${INPUT%.fq.gz}*Met.txt ${INPUT%.fq.gz}*outy "$WD/${INPUT%.fq.gz}_analysis/Bisulfite_alignment"
else
    echo "Final cleanup not engaged, only used with method ALL"
fi


if [ "$METHOD" == "BIN_ONLY" ] ;
# 08/25/19 Removed ALL option, due to the lengthy time this takes to do, WG plots must be engaged by BIN_ONLY
# i.e. nohup sh BSmap_hcDMR.bash /path/to/original.outy /path/to/ref.genome.fa BIN_ONLY &          
then
    echo "Genome binning engaged" ;

    # Separating contexts and making WG methylation plots
    for i in "CG" "CHG" "CHH" ;
    do
        source activate base  # base on my system is python 2.7.15, you'll need to fill in your own

        grep $i ${PRIMARY_INPUT} > ${PRIMARY_INPUT%.outy}.${i}.tmp ;
        
        python $PLOTTER/WG_methylation_plotter.py $PLOTTER/genome_1Mb_binned.txt ${PRIMARY_INPUT%.outy}.${i}.tmp > ${PRIMARY_INPUT%.outy}.${i}.WG_plot ;
        
        if [ -e ${PRIMARY_INPUT%.outy}.${i}.WG_plot ]; then rm ${PRIMARY_INPUT%.outy}.${i}.tmp ; fi        

    done
    # Cleanup
    mkdir -p ${PRIMARY_INPUT%.outy}_analysis/WG_PLOTS ;
    mv ${PRIMARY_INPUT%.outy}*.WG_plot ${PRIMARY_INPUT%.outy}_analysis/WG_PLOTS

else 
    echo "Genome binning skipped, to engage use method BIN_ONLY" ;
fi



### update in future to include: 
# methdiff DMR calling

#python /mnt/ws/home/dsanders/Data/bsmap-2.90/methdiff.py -o ${1%.fq.gz}.CG -l col,${1%.fq.gz} -x CG -d $2 -r 0.4 -b 100 -p 0.01 ${CONTROL_LIB%.fq.gz}.outy ${1%.fq.gz}.outy  ;

#python /mnt/ws/home/dsanders/Data/bsmap-2.90/methdiff.py -o ${1%.fq.gz}.CHG  -l col,${1%.fq.gz} -x CHG -d $2 -r 0.2 -b 100 -p 0.01 ${CONTROL_LIB%.fq.gz}.outy ${1%.fq.gz}.outy ;

#python /mnt/ws/home/dsanders/Data/bsmap-2.90/methdiff.py -o ${1%.fq.gz}.CHH  -l col,${1%.fq.gz} -x CHH -d $2 -r 0.1 -b 100 -p 0.01 ${CONTROL_LIB%.fq.gz}.outy ${1%.fq.gz}.outy ;



### update in future to include: 
# CHH context specific analysis:

#$HOME/anaconda2/bin/python2.7 $HOME/Data/bsmap-2.5/methratio.py -z -u -s /mnt/ws/home/dsanders/Data/bsmap-2.90/samtools/ -o ${1%.fq.gz}.u.2.5.outy -d $2 ${1%.fq.gz}.bam ;

#$HOME/anaconda2/bin/python2.7 $HOME/Data/CHH_context_analysis.py -in ${1%.fq.gz}.u.2.5.outy >> ${1%.fq.gz}.CHH.context.txt ;








# write find out which genes are overlapped by DMRs

#for i in ${1%.fq.gz}.CHH.hyper.DMR ${1%.fq.gz}.CHH.hypo.DMR ;

#do

#/mnt/ws/progs/bin/python2.7 /mnt/ws/home/dsanders/Peak_Overlap_tool/Gene_map_overlap.py /mnt/ws/home/dsanders/Peak_Overlap_tool/TAIR10_genes_only.bed $i -ma 500 > ${i%.DMR}.gene.500.map.txt ;

#done

