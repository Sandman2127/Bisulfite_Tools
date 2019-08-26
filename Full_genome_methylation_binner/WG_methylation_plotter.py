from __future__ import division
import pandas as pd
import sys,os
#from multiprocessing import Pool

#usage python2.7 WG_methylation_plotter.py <binfile see read_bins()> <metratiofile see split metfile>

cwd = os.getcwd()
genels = str(sys.argv[2].strip())
chr1_met = (cwd + "/"+ genels + "_WGplot_chr1.txt")
chr2_met = (cwd + "/"+ genels + "_WGplot_chr2.txt")
chr3_met = (cwd + "/"+ genels + "_WGplot_chr3.txt")
chr4_met = (cwd + "/"+ genels + "_WGplot_chr4.txt")
chr5_met = (cwd + "/"+ genels + "_WGplot_chr5.txt")


def read_bins():

# input file 1: bins
#1	0	500000
#1	500000	1000000

    f = open(sys.argv[1])
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

    with open(sys.argv[2]) as source_file:
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
            next(source_file)
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

def multi_process_wrapper_1(args):
    try:
        return process_chrom(*args)
    except IOError:
        pass

def run_tearDown():
    try:
        os.remove(chr1_met)
        os.remove(chr2_met)
        os.remove(chr3_met)
        os.remove(chr4_met)
        os.remove(chr5_met)
    except OSError:
        pass


if __name__ == "__main__":
    bins = read_bins()

    WG,chr1_df,chr2_df,chr3_df,chr4_df,chr5_df = bins[0],bins[1],bins[2],bins[3],bins[4],bins[5]

    split_metfile()

    #pool = Pool(2)

    #pool.map(multi_process_wrapper_1,[[chr1_df,chr1_met],[chr2_df,chr2_met],[chr3_df,chr3_met],[chr4_df,chr4_met],[chr5_df,chr5_met]])

    process_chrom(chr1_df,chr1_met)

    process_chrom(chr2_df,chr2_met)

    process_chrom(chr3_df,chr3_met)

    process_chrom(chr4_df,chr4_met)

    process_chrom(chr5_df,chr5_met)

    run_tearDown()
