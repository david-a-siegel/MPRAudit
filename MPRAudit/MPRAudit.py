import numpy as np
import os
import scipy.stats as st
import pandas as pd
import argparse

def ratio_function(array1,array2,flag_value=1):
    RNA_counts = np.array(array1)
    DNA_counts = np.array(array2)
    if len(RNA_counts)!=len(DNA_counts):
        raise Exception('RNA_counts and DNA_counts must be array-like of the same length')
    if flag_value==1:
        return float(sum(RNA_counts))/(sum(RNA_counts)+sum(DNA_counts))
    elif flag_value==2:
        return np.log2(float(sum(RNA_counts))/sum(DNA_counts))
    elif flag_value==3:
        return np.log2(float(sum(RNA_counts)+1)/(sum(DNA_counts)+1))
    elif flag_value==4:
        return float(sum(RNA_counts))/sum(DNA_counts)
    elif flag_value==5:
        return (float(sum(RNA_counts)+1))/(sum(DNA_counts)+1)
    else:
        raise Exception('ratio function must have a flag')

def uncertainty(a,b,da,db):
    da = np.sqrt(da) #I'm inputting the variance, should be std dev.
    db = np.sqrt(db) #I'm inputting the variance, shoudl be std dev.
    return 1-a/b,((b*da)**2+(a*db)**2)/b**4


def ordinary_jackknife_variance(list_of_values):
    list_of_values = np.array(list_of_values)
    
    x_i = []
    for i in range(len(list_of_values)):
        list_of_values_i = np.append(list_of_values[:i],list_of_values[i+1:])
        x_i.append(np.mean(list_of_values_i))
    x_i = np.array(x_i)
    jack_var_x = float(len(list_of_values)-1)/len(list_of_values)*sum((x_i-np.mean(x_i))**2)
    return jack_var_x

def ordinary_variance_jackknife_variance(list_of_values):
    list_of_values = np.array(list_of_values)
    
    x_i = []
    for i in range(len(list_of_values)):
        list_of_values_i = np.append(list_of_values[:i],list_of_values[i+1:])
        x_i.append(np.var(list_of_values_i))
    x_i = np.array(x_i)
    jack_var_x = float(len(list_of_values)-1)/len(list_of_values)*sum((x_i-np.mean(x_i))**2)
    return jack_var_x
    
    
def ordinary_jackknife_variance_Ngroups(list_of_values,group_indicators):
    list_of_values = np.array(list_of_values)
    group_indicators = np.array(group_indicators)
    N = len(list_of_values)
    
    x_i = []
    for group in np.unique(group_indicators):
        group_values = list_of_values[group_indicators==group]
        notgroup_values = list_of_values[group_indicators!=group]
        for i in range(len(group_values)):
            group_values_i = np.append(group_values[:i],group_values[i+1:])
            x_i.append(np.mean(np.append(group_values_i,notgroup_values)))
    
    x_i = np.array(x_i)
    jack_var_x = float(N-1)/N*sum((x_i-np.mean(x_i))**2)
    return jack_var_x


def ordinary_variance_jackknife_variance_Ngroups(list_of_values,group_indicators):
    list_of_values = np.array(list_of_values)
    group_indicators = np.array(group_indicators)
    N = len(list_of_values)

    x_i = []
    for group in np.unique(group_indicators):
        group_values = list_of_values[group_indicators==group]
        for i in range(len(group_values)):
            group_values_i = np.append(group_values[:i],group_values[i+1:])
            running_sum = 0.
            for group_prime in np.unique(group_indicators):
                if group_prime == group:
                    running_sum += len(group_values_i)*np.var(group_values_i)
                else:
                    running_sum += len(list_of_values[group_indicators==group_prime])*np.var(list_of_values[group_indicators==group_prime])
            x_i.append(running_sum/(N-1))
        
    
    x_i = np.array(x_i)
    jack_var_x = float(N-1)/N*sum((x_i-np.mean(x_i))**2)
    return jack_var_x



def deleteDjackknife_variance_T0(RNA_counts, DNA_counts, num_trials = 100, exp_pow = 3./5, ratio_flag = 1):
    if len(RNA_counts)!=len(DNA_counts):
        print("Error: RNA counts has different length from DNA counts")
    RNA_counts = np.array(RNA_counts)
    DNA_counts = np.array(DNA_counts)

    d = int(len(RNA_counts)**(exp_pow))
    d_kept = len(RNA_counts) - d #This is n-d
    
    x_i = []
    for i in range(num_trials):
        kept_clones = random.sample(range(len(RNA_counts)),d_kept)
        RNA_counts_i = RNA_counts[kept_clones]
        DNA_counts_i = DNA_counts[kept_clones]
        x_i.append(ratio_function(RNA_counts_i,DNA_counts_i,ratio_flag))
    x_i = np.array(x_i)
    jack_var_x = float(d_kept)/d/num_trials*sum((x_i-np.mean(x_i))**2)
    return jack_var_x

    
    
def MPRAudit_Pairs(RNA_counts1, DNA_counts1, RNA_counts2, DNA_counts2, sequence_indicators, numtrials, jackpow, ratio_flag) #sequence_indicator groups clones into otherwise identical sequences
    if len(RNA_counts1)!=len(DNA_counts1):
        raise Exception('RNA_counts and DNA_counts must be array-like of the same length')
    if len(RNA_counts1)!=len(RNA_counts2):
        raise Exception('RNA_counts1 and RNA_counts2 must be array-like of the same length')
    if len(RNA_counts1)!=len(DNA_counts2):
        raise Exception('RNA_counts and DNA_counts must be array-like of the same length')        
    if len(sequence_indicators)!=len(RNA_counts):
        raise Exception('sequence_indicators must be array-like of the same length as the counts')

    RNA_counts1 = np.array(RNA_counts1)
    DNA_counts1 = np.array(DNA_counts1)
    RNA_counts2 = np.array(RNA_counts2)
    DNA_counts2 = np.array(DNA_counts2)
    sequence_indicators = np.array(sequence_indicators)
    
    
    jackknife_variance_list = []
    RD_list = []
    
    for sequence in np.unique(sequence_indicators):
        RNA_sequence_counts1 = RNA_counts1[sequence_indicators==sequence]
        DNA_sequence_counts1 = DNA_counts1[sequence_indicators==sequence]    
        RNA_sequence_counts2 = RNA_counts2[sequence_indicators==sequence]
        DNA_sequence_counts2 = DNA_counts2[sequence_indicators==sequence]    
        jackknife_variance_list.append(deleteDjackknife_variance_T0(RNA_sequence_counts1, DNA_sequence_counts1,numtrials,jackpow,ratio_flag)+deleteDjackknife_variance_T0(RNA_sequence_counts2, DNA_sequence_counts2,numtrials,jackpow,ratio_flag))
        RD_list.append(ratio_function(RNA_sequence_counts1,DNA_sequence_counts1,ratio_flag)-ratio_function(RNA_sequence_counts2,DNA_sequence_counts2,ratio_flag))
        
    tech_variance = np.mean(jackknife_variance_list)
    total_variance = np.var(RD_list)
    delta_tech = ordinary_jackknife_variance(jackknife_variance_list)
    delta_var = ordinary_variance_jackknife_variance(RD_list))
    b2_mean, b2_var = uncertainty(tech_variance,total_variance,delta_tech,delta_var)
    return b2_mean, np.sqrt(b2_var)
    
    
def MPRAudit_Groups(RNA_counts, DNA_counts, sequence_indicators, numtrials, jackpow, ratiofunction, group_indicators=None) #sequence_indicator groups clones into otherwise identical sequences
    if len(RNA_counts)!=len(DNA_counts):
        raise Exception('RNA_counts and DNA_counts must be array-like of the same length')
    if len(sequence_indicators)!=len(RNA_counts):
        raise Exception('sequence_indicators must be array-like of the same length as the counts')
    if group_indicators is not None:
        if len(group_indicators)!=len(RNA_counts):
            raise Exception('sequence_indicators must be array-like of the same length as the counts')    

    RNA_counts = np.array(RNA_counts)
    DNA_counts = np.array(DNA_counts)
    sequence_indicators = np.array(sequence_indicators)
    group_indicators = np.array(group_indicators)
    #Check lengths?
    
    tech_variance = 0
    total_variance = 0
    for group in np.unique(group_indicators):
        jackknife_variance_list = []
        RD_list = []
        for sequence in np.unique(sequence_indicators):
            RNA_sequence_counts = RNA_counts[(sequence_indicators==sequence)&(group_indicators==group)]
            DNA_sequence_counts = DNA_counts[(sequence_indicators==sequence)&(group_indicators==group)]
            jackknife_variance_list.append(deleteDjackknife_variance_T0(RNA_sequence_counts, DNA_sequence_counts,numtrials,jackpow,ratio_flag))
            RD_list.append(ratio_function(RNA_sequence_counts,DNA_sequence_counts,ratio_flag))
        tech_variance += np.mean(jackknife_variance_list)*len(jackknife_variance_list)
        total_variance += np.var(RD_list)*len(RD_list)
    tech_variance = tech_variance/len(RNA_counts)
    total_variance = total_variance/len(RNA_counts)
    delta_tech = ordinary_jackknife_variance_Ngroups(jackknife_variance_list,group_indicators)
    delta_var = ordinary_variance_jackknife_variance_Ngroups(RD_list,group_indicators)
    b2_mean, b2_var = uncertainty(tech_variance,total_variance,delta_tech,delta_var)
    return b2_mean, np.sqrt(b2_var)
    

def MPRAudit_function(infile,outfile,ratiofunction,paired,sepstr,timepoints,numtrials,jackpow):

    data_DF = pd.read_csv(infile,sep=sepstr,header=None)
    if timepoints==1:
        if paired==False:
            if data_DF.shape[1]==3:
                b2_mean,b2_std = MPRAudit_Groups(data_DF.iloc[:,0],data_DF.iloc[:,1],data_DF.iloc[:,2],numtrials=numtrials,jackpow=jackpow,ratiofunction=ratiofunction) #No groups, just single sequences
            elif data_DF.shape[1]==4:    
                b2_mean,b2_std = MPRAudit_Groups(data_DF.iloc[:,0],data_DF.iloc[:,1],data_DF.iloc[:,2],data_DF.iloc[:,3],numtrials=numtrials,jackpow=jackpow,ratiofunction=ratiofunction) #Groups
            else:
                raise Exception('Input file should have 3 or 4 columns for this set of parameters')
        elif paired==True: #Not single sequences, difference between sequences
            if data_DF.shape[1]==5:
                b2_mean,b2_std = MPRAudit_Pairs(data_DF.iloc[:,0],data_DF.iloc[:,1],data_DF.iloc[:,2],data_DF.iloc[:,3],data_DF.iloc[:,4],numtrials=numtrials,jackpow=jackpow,ratiofunction=ratiofunction)
            else:
                raise Exception('Input file should have 5 columns for this set of parameters')

    return b2_mean, b2_std    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-infile", "--infile", type=str, help="Input filename")
    parser.add_argument("-outfile", "--outfile",type=str, default="outfile",help="Output filename, outfile by default")
    parser.add_argument("-ratio","--ratiofunction",type=int,default=1,help="Which ratio function to use, see help file for details")
    parser.add_argument("-paired","--paired",type=bool,default=False,help="Paired sequences, default = False")
    parser.add_argument("-sep","--sepstr",type=str,default="\t",help="CSV or txt delimiter, default = \t"
    parser.add_argument("-timepoints","--timepoints",type=int,default=1,help="T0 or T4/(T4+T0)")
    parser.add_argument("-numtrials","--numtrials",type=int,default=100,help="Number of jackknife trials")
    parser.add_argument("-jackpow","--jackpow",type=float,default=3./5,help = "Jackknife power exponent")
    

    b2_mean, b2_std = MPRAudit_function(args.infile,args.outfile,args.ratiofunction,args.paired,args.sepstr,args.timepoints,args.numtrials,args.jackpow)
    
    
    pd.DataFrame(data={"b2_mean":b2_mean,"b2_std":b2_std}).to_csv(outfile)
                            

if __name__ == '__main__':
     main()

