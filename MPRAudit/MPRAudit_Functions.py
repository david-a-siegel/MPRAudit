import numpy as np
import os
import scipy.stats as st
import pandas as pd
import random

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


def generate_random_seq(length_of_seq,has_TATACAG):
    seq_vals = np.random.randint(0,4,size=length_of_seq)
    seq_string = ''
    for i in range(length_of_seq):
        if seq_vals[i]==0:
            seq_string+='A'
        elif seq_vals[i]==1:
            seq_string+='C'
        elif seq_vals[i]==2:
            seq_string+='T'
        elif seq_vals[i]==3:
            seq_string+='G'
        else:
            print('error in string generation')
    if has_TATACAG:
        start_position = np.random.randint(0,length_of_seq-len('TATACAG'))
        seq_string = seq_string[:start_position]+'TATACAG'+seq_string[start_position+len('TATACAG'):]
    return seq_string
    

def return_counts(num_clones=100,dna_per_clone=200,A=3,s=0.1,sigma_mu_ratio=30.,K=1.):
    DNA_true = st.poisson.rvs(dna_per_clone,size=num_clones)
    DNA_true[DNA_true==0]=1
    RNA_true = np.round(A*DNA_true + st.norm.rvs(scale=s*np.sqrt(DNA_true),size=num_clones))
    RNA_true[RNA_true<0]=0 #Make sure you don't have negative number of transcripts

    KR = 1.25*K
    muRNA = RNA_true*KR
    sigma2R = muRNA*sigma_mu_ratio+1e-6
    rRNA = muRNA**2/(sigma2R-muRNA)
    pRNA = rRNA*1.0/(rRNA+muRNA)
    RNA_sampled = []
    for i in range(len(rRNA)):
        if rRNA[i]==0:
            RNA_sampled.append(0)
        else:
            RNA_sampled.append(st.nbinom.rvs(rRNA[i],pRNA[i]))
            
    KD = 1.0*K
    muDNA = DNA_true*KD
    sigma2D = muDNA*sigma_mu_ratio
    rDNA = muDNA**2/(sigma2D-muDNA)
    pDNA = rDNA*1.0/(rDNA+muDNA)
    DNA_sampled = []
    for i in range(len(DNA_true)):
        if rDNA[i]==0:
            DNA_sampled.append(0)
        else:
            DNA_sampled.append(st.nbinom.rvs(rDNA[i],pDNA[i]))
    
    return RNA_sampled,DNA_sampled, KR, KD
    
    
def return_counts_T4T0(num_clones=100,dna_per_clone=200,A=3,s=1.,sigma_mu_ratio=30.,T4T0_ratio=0.5,K=1):
    DNA_true_T0 = st.poisson.rvs(dna_per_clone,size=num_clones)
    DNA_true_T0[DNA_true_T0==0]=1

    DNA_true_T4 = st.poisson.rvs(dna_per_clone,size=num_clones)
    DNA_true_T4[DNA_true_T4==0]=1
    
    RNA_true_T0 = np.round(A*DNA_true_T0 + st.norm.rvs(scale=s*np.sqrt(DNA_true_T0),size=num_clones))
    RNA_true_T0[RNA_true_T0<0]=0 #Make sure RNA counts>=0

    RNA_true_T4 = np.round(A*DNA_true_T4*T4T0_ratio+st.norm.rvs(scale=s*np.sqrt(DNA_true_T4),size=num_clones))
    RNA_true_T4[RNA_true_T4<0]=0 #Make sure RNA counts>=0

    #Make them slightly different to avoid accidental cancellation by symmetry    
    KR0 = 1.3*K
    KD0 = 1.0*K
    KR4 = 2.2*K
    KD4 = 2.4*K


    
    muRNA_T0 = RNA_true_T0*KR0
    sigma2R_T0 = muRNA_T0*sigma_mu_ratio+1e-6
    rRNA_T0 = muRNA_T0**2/(sigma2R_T0-muRNA_T0)
    pRNA_T0 = rRNA_T0*1.0/(rRNA_T0+muRNA_T0)

    muRNA_T4 = RNA_true_T4*KR4
    sigma2R_T4 = muRNA_T4*sigma_mu_ratio+1e-6
    rRNA_T4 = muRNA_T4**2/(sigma2R_T4-muRNA_T4)
    pRNA_T4 = rRNA_T4*1.0/(rRNA_T4+muRNA_T4)

    
    RNA_sampled_T0 = []
    for i in range(len(rRNA_T0)):
        if rRNA_T0[i]==0:
            RNA_sampled_T0.append(0)
        else:
            RNA_sampled_T0.append(st.nbinom.rvs(rRNA_T0[i],pRNA_T0[i]))
    
    RNA_sampled_T4 = []
    for i in range(len(rRNA_T4)):
        if rRNA_T4[i]==0:
            RNA_sampled_T4.append(0)
        else:
            RNA_sampled_T4.append(st.nbinom.rvs(rRNA_T4[i],pRNA_T4[i]))
            

    muDNA_T0 = DNA_true_T0*KD0
    sigma2D_T0 = muDNA_T0*sigma_mu_ratio
    rDNA_T0 = muDNA_T0**2/(sigma2D_T0-muDNA_T0)
    pDNA_T0 = rDNA_T0*1.0/(rDNA_T0+muDNA_T0)
    
    muDNA_T4 = DNA_true_T4*KD4
    sigma2D_T4 = muDNA_T4*sigma_mu_ratio
    rDNA_T4 = muDNA_T4**2/(sigma2D_T4-muDNA_T4)
    pDNA_T4 = rDNA_T4*1.0/(rDNA_T4+muDNA_T4)

    
    DNA_sampled_T0 = []
    for i in range(len(DNA_true_T0)):
        if rDNA_T0[i]==0:
            DNA_sampled_T0.append(0)
        else:
            DNA_sampled_T0.append(st.nbinom.rvs(rDNA_T0[i],pDNA_T0[i]))

    DNA_sampled_T4 = []
    for i in range(len(DNA_true_T4)):
        if rDNA_T4[i]==0:
            DNA_sampled_T4.append(0)
        else:
            DNA_sampled_T4.append(st.nbinom.rvs(rDNA_T4[i],pDNA_T4[i]))
            
            
    return RNA_sampled_T0,DNA_sampled_T0,RNA_sampled_T4,DNA_sampled_T4, KR0, KD0, KR4, KD4

    
    

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



def deleteDjackknife_variance_T0(RNA_counts, DNA_counts, num_trials = 100, exp_pow = 3./5, ratiofunction = 1):
    if len(RNA_counts)!=len(DNA_counts):
        raise Exception("Error: RNA counts has different length from DNA counts")
    RNA_counts = np.array(RNA_counts)
    DNA_counts = np.array(DNA_counts)

    d = int(len(RNA_counts)**(exp_pow))
    d_kept = len(RNA_counts) - d #This is n-d
    
    x_i = []
    for i in range(num_trials):
        kept_clones = random.sample(range(len(RNA_counts)),d_kept)
        RNA_counts_i = RNA_counts[kept_clones]
        DNA_counts_i = DNA_counts[kept_clones]
        x_i.append(ratio_function(RNA_counts_i,DNA_counts_i,ratiofunction))
    x_i = np.array(x_i)
    jack_var_x = float(d_kept)/d/num_trials*sum((x_i-np.mean(x_i))**2)
    return jack_var_x

    
    
def MPRAudit_Pairs(RNA_counts1, DNA_counts1, RNA_counts2, DNA_counts2, sequence_indicators1, sequence_indicators2, numtrials, jackpow, ratiofunction): #sequence_indicator groups clones into otherwise identical sequences
    if len(RNA_counts1)!=len(DNA_counts1):
        raise Exception('RNA_counts and DNA_counts must be array-like of the same length')
    if len(RNA_counts2)!=len(DNA_counts2):
        raise Exception('RNA_counts and DNA_counts must be array-like of the same length')        
    if len(sequence_indicators1)!=len(RNA_counts1):
        raise Exception('sequence_indicators must be array-like of the same length as the counts')
    if len(sequence_indicators2)!=len(RNA_counts2):
        raise Exception('sequence_indicators must be array-like of the same length as the counts')
    

    RNA_counts1 = np.array(RNA_counts1)
    DNA_counts1 = np.array(DNA_counts1)
    RNA_counts2 = np.array(RNA_counts2)
    DNA_counts2 = np.array(DNA_counts2)
    sequence_indicators1 = np.array(sequence_indicators1)
    sequence_indicators2 = np.array(sequence_indicators2)
    RNA_counts1 = RNA_counts1[~pd.isnull(RNA_counts1)]
    RNA_counts2 = RNA_counts2[~pd.isnull(RNA_counts2)]
    DNA_counts1 = DNA_counts1[~pd.isnull(DNA_counts1)]
    DNA_counts2 = DNA_counts2[~pd.isnull(DNA_counts2)]
    sequence_indicators1 = sequence_indicators1[~pd.isnull(sequence_indicators1)]
    sequence_indicators2 = sequence_indicators2[~pd.isnull(sequence_indicators2)]    
    
    jackknife_variance_list = []
    RD_list = []
    
    if set(sequence_indicators1)!=set(sequence_indicators2):
        raise Exception('Sequence indicators must be paired in MPRAudit Paired')
        
    for sequence in np.unique(sequence_indicators1):
        RNA_sequence_counts1 = RNA_counts1[sequence_indicators1==sequence]
        DNA_sequence_counts1 = DNA_counts1[sequence_indicators1==sequence]    
        RNA_sequence_counts2 = RNA_counts2[sequence_indicators2==sequence]
        DNA_sequence_counts2 = DNA_counts2[sequence_indicators2==sequence]    
        jackknife_variance_list.append(deleteDjackknife_variance_T0(RNA_sequence_counts1, DNA_sequence_counts1,numtrials,jackpow,ratiofunction)+deleteDjackknife_variance_T0(RNA_sequence_counts2, DNA_sequence_counts2,numtrials,jackpow,ratiofunction))
        RD_list.append(ratio_function(RNA_sequence_counts1,DNA_sequence_counts1,ratiofunction)-ratio_function(RNA_sequence_counts2,DNA_sequence_counts2,ratiofunction))
        
    tech_variance = np.mean(jackknife_variance_list)
    total_variance = np.var(RD_list)
    delta_tech = ordinary_jackknife_variance(jackknife_variance_list)
    delta_var = ordinary_variance_jackknife_variance(RD_list)
    b2_mean, b2_var = uncertainty(tech_variance,total_variance,delta_tech,delta_var)
    return b2_mean, np.sqrt(b2_var)
    
    
def MPRAudit_Groups(RNA_counts, DNA_counts, sequence_indicators, numtrials=100, jackpow=3./5, ratiofunction=1, group_indicators=None): #sequence_indicator groups clones into otherwise identical sequences
    if len(RNA_counts)!=len(DNA_counts):
        raise Exception('RNA_counts and DNA_counts must be array-like of the same length')
    if len(sequence_indicators)!=len(RNA_counts):
        raise Exception('sequence_indicators must be array-like of the same length as the counts')
    if group_indicators is not None:
        if len(group_indicators)!=len(RNA_counts):
            raise Exception('sequence_indicators must be array-like of the same length as the counts')    
    else: #If no group indicators then no groups -- just make them all one.
        group_indicators = np.ones(len(RNA_counts))
    RNA_counts = np.array(RNA_counts)
    DNA_counts = np.array(DNA_counts)
    sequence_indicators = np.array(sequence_indicators)
    group_indicators = np.array(group_indicators)
    
    tech_variance = 0
    total_variance = 0
    group_ordered_list = []
    
    for group in np.unique(group_indicators):
        jackknife_variance_list = []
        RD_list = []
        for sequence in np.unique(sequence_indicators):
            RNA_sequence_counts = RNA_counts[(sequence_indicators==sequence)&(group_indicators==group)]
            DNA_sequence_counts = DNA_counts[(sequence_indicators==sequence)&(group_indicators==group)]
            jackknife_variance_list.append(deleteDjackknife_variance_T0(RNA_sequence_counts, DNA_sequence_counts,numtrials,jackpow,ratiofunction))
            RD_list.append(ratio_function(RNA_sequence_counts,DNA_sequence_counts,ratiofunction))
            if len(RD_list)>0:
                group_ordered_list.append(group) #I need the group indicators in the right order for later, one for each sequence now instead of one for each clone
        tech_variance += np.mean(jackknife_variance_list)*len(jackknife_variance_list)
        total_variance += np.var(RD_list)*len(RD_list)

    tech_variance = tech_variance/len(group_ordered_list) #Divide by total number of sequences
    total_variance = total_variance/len(group_ordered_list)
    delta_tech = ordinary_jackknife_variance_Ngroups(jackknife_variance_list,np.array(group_ordered_list))
    delta_var = ordinary_variance_jackknife_variance_Ngroups(RD_list,np.array(group_ordered_list))

    b2_mean, b2_var = uncertainty(tech_variance,total_variance,delta_tech,delta_var)
    return b2_mean, np.sqrt(b2_var)
    

def MPRAudit_function(data_DF,ratiofunction=1,paired=False,timepoints=1,numtrials=100,jackpow=3./5):

    if timepoints==1:
        if paired==False:
            if data_DF.shape[1]==3:
                b2_mean,b2_std = MPRAudit_Groups(data_DF.iloc[:,0],data_DF.iloc[:,1],data_DF.iloc[:,2],numtrials=numtrials,jackpow=jackpow,ratiofunction=ratiofunction) #No groups, just single sequences
            elif data_DF.shape[1]==4:    
                b2_mean,b2_std = MPRAudit_Groups(data_DF.iloc[:,0],data_DF.iloc[:,1],data_DF.iloc[:,2],group_indicators = data_DF.iloc[:,3],numtrials=numtrials,jackpow=jackpow,ratiofunction=ratiofunction) #Groups
            else:
                raise Exception('Input file should have 3 or 4 columns for this set of parameters')
        elif paired==True: #Not single sequences, difference between sequences
            if data_DF.shape[1]==6:
                b2_mean,b2_std = MPRAudit_Pairs(data_DF.iloc[:,0],data_DF.iloc[:,1],data_DF.iloc[:,2],data_DF.iloc[:,3],data_DF.iloc[:,4],data_DF.iloc[:,5],numtrials=numtrials,jackpow=jackpow,ratiofunction=ratiofunction)
            else:
                raise Exception('Input file should have 5 columns for this set of parameters')

    return b2_mean, b2_std    
