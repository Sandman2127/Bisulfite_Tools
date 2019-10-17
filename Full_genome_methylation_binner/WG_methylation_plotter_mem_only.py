from __future__ import division
import pandas as pd
import sys,os
import argparse

parser = argparse.ArgumentParser(prog='Whole Genome DNA Methylation Plotter')
parser.add_argument('-input', dest='METFILE',action="store", type=argparse.FileType('r'),help="tab delimited output (.outy) file BSmap Methratio.py")
parser.add_argument('-binfile', dest='BINFILE',action="store", type=argparse.FileType('r'),help="tab delimited bin file i.e. :    1  0    500000")

args = parser.parse_args()
_METFILE=args.METFILE.name
_BINFILE=args.BINFILE.name
cwd = os.getcwd()


def read_bins():

# input file 1: bins
#1	0	500000
#1	500000	1000000

    f = open(_BINFILE)
    df = pd.read_table(f, sep='\t',names=('chromosome','start','end'))
    df_chr1 = df[df['chromosome'] == 1]
    df_chr2 = df[df['chromosome'] == 2]
    df_chr3 = df[df['chromosome'] == 3]
    df_chr4 = df[df['chromosome'] == 4]
    df_chr5 = df[df['chromosome'] == 5]
    list_o_dfs = [df,df_chr1,df_chr2,df_chr3,df_chr4,df_chr5]
    return list_o_dfs
    f.close()

def split_metfile():

#file 2: stuff to bin
#chr1	8	+	CHH	0.000	1.00	0	1	3	3	0.000	0.793
#chr1	9	+	CHH	1.000	1.00	1	1	4	4	0.207	1.000
#chr1	10	+	CHH	1.000	2.00	2	2	4	4	0.342	1.000

    with open(_METFILE) as source_file:
        next(source_file) 
        for line in source_file:  # begin iteration over the input file line by line
            line = line.strip() # strip of whitespace
            cyt = line.split()   # this splits all the words into columns
            cyt_chr = cyt[0].strip()
            if cyt_chr == "chr1":
                with open(chr1_met, 'a+') as f:   # a allows you to append to that file
                    f.write(str(cyt[0][-1:]) + "\t" + str(cyt[1]) + "\t" + str(cyt[2]) + "\t" + str(cyt[3]) + "\t" + str(cyt[4])+"\n")
            if cyt_chr == "chr2":
                with open(chr2_met, 'a+') as f:   # a allows you to append to that file
                    f.write(str(cyt[0][-1:]) + "\t" + str(cyt[1]) + "\t" + str(cyt[2]) + "\t" + str(cyt[3]) + "\t" + str(cyt[4])+"\n")
            if cyt_chr == "chr3":
                with open(chr3_met, 'a+') as f:   # a allows you to append to that file
                    f.write(str(cyt[0][-1:]) + "\t" + str(cyt[1]) + "\t" + str(cyt[2]) + "\t" + str(cyt[3]) + "\t" + str(cyt[4])+"\n")
            if cyt_chr == "chr4":
                with open(chr4_met, 'a+') as f:   # a allows you to append to that file
                    f.write(str(cyt[0][-1:]) + "\t" + str(cyt[1]) + "\t" + str(cyt[2]) + "\t" + str(cyt[3]) + "\t" + str(cyt[4])+"\n")
            if cyt_chr == "chr5":
                with open(chr5_met, 'a+') as f:   # a allows you to append to that file
                    f.write(str(cyt[0][-1:]) + "\t" + str(cyt[1]) + "\t" + str(cyt[2]) + "\t" + str(cyt[3]) + "\t" + str(cyt[4])+"\n")

# temp write file format:
#1	2	+	CHH	0.000
#1	3	+	CHH	0.000
#1	8	+	CHH	0.000

def process_chrom(bins,chr_file):
    index = 0
    for row in bins.itertuples():
        list1 = []
        index += 1
        start = int(row[2])
        end = int(row[3])
        binratio = 0
        bintotal = 1  # makes sure each bin has at least 1 C
        with open(chr_file) as source_file:
            next(source_file)  #probably skips the first line 
            for line in source_file:  # begin iteration over the input file line by line
                line = line.strip()   # strips the line of white space
                cyt = line.split()   # this splits all the words into columns
                cyt_pos,ratio = int(cyt[1]),float(cyt[4])
                if cyt_pos >= start and cyt_pos <= end:
                    binratio += ratio
                    bintotal += 1
                elif cyt_pos > end:
                    break
                else:
                    pass
        print index,"\t",100*(binratio/bintotal),"\t",start,"\t",end

def return_met_conversion_rate():
    with open(_METFILE) as outyfile:
        inputdf = pd.read_table(outyfile,sep='\t')
        
                #shape of full, few additional columns not used
                #     chr  pos strand context  ratio  eff_CT_count  C_count  CT_count  \
                # 0  chr1   15      +     CHH    0.0           1.0        0         1   
                # 1  chr1   16      +     CHH    0.0           1.0        0         1   
                # 2  chr1   17      +     CHH    1.0           1.0        1         1   
                # 3  chr1   22      +     CHH    0.0           1.0        0         1   
                # 4  chr1   23      +     CHH    0.0           1.0        0         1

        chrC = inputdf[inputdf['chr'] == 'chrC']
        listA=chrC.sum(axis=0,skipna=True).tolist()
        C_count=float(listA[6])
        CT_count=float(listA[7])
        Unconv=100*(C_count/CT_count)
        Conv_efficiency=100-Unconv
        print "Conversion efficiency equals: ",Conv_efficiency," normalized results will subtract ",Unconv,"from the raw total"
        return float(Unconv)



def split_metfile_df(context):
    with open(_METFILE) as outyfile:
        full = pd.read_table(outyfile,sep='\t')
        
                #shape of full, few additional columns not used
                #     chr  pos strand context  ratio  eff_CT_count  C_count  CT_count  \
                # 0  chr1   15      +     CHH    0.0           1.0        0         1   
                # 1  chr1   16      +     CHH    0.0           1.0        0         1   
                # 2  chr1   17      +     CHH    1.0           1.0        1         1   
                # 3  chr1   22      +     CHH    0.0           1.0        0         1   
                # 4  chr1   23      +     CHH    0.0           1.0        0         1

        chr1 = full[(full['chr'] == 'chr1') & (full['context'] == str(context))]
        chr2 = full[(full['chr'] == 'chr2') & (full['context'] == str(context))]
        chr3 = full[(full['chr'] == 'chr3') & (full['context'] == str(context))]
        chr4 = full[(full['chr'] == 'chr4') & (full['context'] == str(context))]
        chr5 = full[(full['chr'] == 'chr5') & (full['context'] == str(context))]
        next_list_o_dfs = [chr1,chr2,chr3,chr4,chr5]
        return next_list_o_dfs


def process_chrom_df(bins,met_df,corr_factor,writefile):
    #TODO process methylation into bins 
    index = 0
    for row in bins.itertuples():
        list1 = []
        index += 1
        start = int(row[2])
        end = int(row[3])
        binratio = float(0)
        bintotal = 1 
        for metCs in met_df.itertuples():
            cyt_pos = int(metCs[2])
            ratio = float(metCs[5])
            if cyt_pos >= start and cyt_pos <= end:
                binratio += ratio
                bintotal += 1
            elif cyt_pos > end:
                break
            else:
                pass
        # protects for situations where this value may go negative due to incorrect corr_factor or other methylation differences
        efficiency_corrected_methylation=(100*(binratio/bintotal))-corr_factor
        if efficiency_corrected_methylation < 0: 
            output=str(index) + "\t" + str(100*(binratio/bintotal)) + "\t" + str(0) + "\t" + str(start) + "\t" + str(end) + "\n"
            writefile.write(output)
        else:
            output=str(index) + "\t" + str(100*(binratio/bintotal)) + "\t" + str(efficiency_corrected_methylation) + "\t" + str(start) + "\t" + str(end) + "\n"
            writefile.write(output)


def multi_process_wrapper_1(args):
    try:
        return process_chrom_df(*args)
        #return process_chrom(*args) #old non-df version
    except IOError:
        pass

if __name__ == "__main__":

    #TODO calc conversion efficiency to use for normalization:
    corr_factor=return_met_conversion_rate()

    #TODO readin bins:
    print "Reading in genome bins..."
    bins = read_bins()
    WG,chr1_bins,chr2_bins,chr3_bins,chr4_bins,chr5_bins = bins[0],bins[1],bins[2],bins[3],bins[4],bins[5]
    print "Bin read COMPLETE"

    #TODO process chromsomes by context:
    context_list = ['CG','CHG','CHH']
    for context in context_list:
        print "Starting analysis of: ",context," context"

        #TODO readin metfile:
        print "Spliting metfile by chromosome and ",context," context"
        chrlist=split_metfile_df(context)
        met1,met2,met3,met4,met5=chrlist[0],chrlist[1],chrlist[2],chrlist[3],chrlist[4]
        
        #TODO open a file to push WG_plot data to:
        fileoutname = cwd + "/" + str(_METFILE.replace(".outy",'')) + "." + str(context) + ".WGplot.txt"
        with open(fileoutname,"w+") as writefile:

            #TODO write header
            header=str("Chrom_Index") + "\t" + str("Uncorrected_Methylation") + "\t" + str("Conv_Eff_Corrected_Methylation") + "\t" + str("Bin_Start") + "\t" + str("Bin_End") + "\n"
            writefile.write(header)

            #TODO process in this context:
            print "Writing binned whole genome methylation data into",fileoutname 
            print "Processing chr1"
            process_chrom_df(chr1_bins,met1,corr_factor,writefile)
            print "Processing chr2"
            process_chrom_df(chr2_bins,met2,corr_factor,writefile)
            print "Processing chr3"
            process_chrom_df(chr3_bins,met3,corr_factor,writefile)
            print "Processing chr4"
            process_chrom_df(chr4_bins,met4,corr_factor,writefile)
            print "Processing chr5"
            process_chrom_df(chr5_bins,met5,corr_factor,writefile)
            print "Data processing complete for: ",fileoutname



