import numpy as np
import os
import scipy.stats as st
import pandas as pd
import argparse
import MPRAudit_Functions

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-infile", "--infile", type=str, help="Input filename")
    parser.add_argument("-outfile", "--outfile",type=str, default=None,help="Output filename, if None then prints to standard output")
    parser.add_argument("-ratio","--ratiofunction",type=int,default=1,help="Which ratio function to use, see help file for details")
    parser.add_argument("-paired","--paired",type=bool,default=False,help="Paired sequences, default = False")
    parser.add_argument("-sep","--sepstr",type=str,default=",",help="CSV or txt delimiter, default = ','")
    parser.add_argument("-timepoints","--timepoints",type=int,default=1,help="T0 or T4/(T4+T0)")
    parser.add_argument("-numtrials","--numtrials",type=int,default=100,help="Number of jackknife trials")
    parser.add_argument("-jackpow","--jackpow",type=float,default=3./5,help = "Jackknife power exponent")
    parser.add_argument("-CRISPR_log_flag","--CRISPR_log_flag",type=bool,default=None,help = "Should be None unless RNA-data only and no DNA data; whether or not data should be log-normalized by MPRAudit")
    args = parser.parse_args()
    
    data_DF = pd.read_csv(args.infile,sep=args.sepstr,header=None) #No header!
    b2_mean, b2_std = MPRAudit_Functions.MPRAudit_function(data_DF,args.ratiofunction,args.paired,args.timepoints,args.numtrials,args.jackpow,args.CRISPR_log_flag)
    
    if args.outfile is None:
        print("b2_mean: "+repr(b2_mean))
        print("b2_std: "+repr(b2_std))
    else:
        pd.DataFrame(data={"b2_mean":b2_mean,"b2_std":b2_std}).to_csv(outfile)
                            

if __name__ == '__main__':
     main()

