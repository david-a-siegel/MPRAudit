import numpy as np
import os
import scipy.stats as st
import pandas as pd
import MPRAudit_Functions
from time import process_time


if __name__ == '__main__':
    starttime = process_time()
    #Let's make sets of pairs.  1 vs 2, where half of 2 is different from the other half.
    RNA_data1 = np.array([])
    DNA_data1 = np.array([])
    RNA_data2 = np.array([])
    DNA_data2 = np.array([])

    sequence_indicators1 = np.array([])
    for sequence_number in range(500):
        RNA_sampled, DNA_sampled, sum_RNA, sum_DNA = MPRAudit_Functions.return_counts(num_clones=20,dna_per_clone=200,A=abs(np.random.normal(loc=3,scale=.1))+0.1,s=0.1,sigma_mu_ratio=30.,K=10.)
        RNA_sampled = np.array(RNA_sampled)/sum_RNA #Save as numpy array and count-normalize (make sure total counts (K) are held constant!)
        DNA_sampled = np.array(DNA_sampled)/sum_DNA
        
        RNA_data1 = np.append(RNA_data1,RNA_sampled)
        DNA_data1 = np.append(DNA_data1,DNA_sampled)
        sequence_indicators1 = np.append(sequence_indicators1,sequence_number*np.ones(len(RNA_sampled)))
    all_zeros = (RNA_data1==0)&(DNA_data1==0)
    RNA_data1 = RNA_data1[~all_zeros]
    DNA_data1 = DNA_data1[~all_zeros]
    sequence_indicators1 = sequence_indicators1[~all_zeros]

    sequence_indicators2 = np.array([])
    for sequence_number in range(250):
        RNA_sampled, DNA_sampled, sum_RNA, sum_DNA = MPRAudit_Functions.return_counts(num_clones=20,dna_per_clone=200,A=abs(np.random.normal(loc=3,scale=1))+0.1+1,s=0.1,sigma_mu_ratio=30.,K=1.)
        RNA_sampled = np.array(RNA_sampled)/sum_RNA #Save as numpy array and count-normalize (make sure total counts (K) are held constant!)
        DNA_sampled = np.array(DNA_sampled)/sum_DNA
        
        RNA_data2 = np.append(RNA_data2,RNA_sampled)
        DNA_data2 = np.append(DNA_data2,DNA_sampled)
        sequence_indicators2 = np.append(sequence_indicators2,sequence_number*np.ones(len(RNA_sampled)))

    for sequence_number in range(250,500):
        RNA_sampled, DNA_sampled, sum_RNA, sum_DNA = MPRAudit_Functions.return_counts(num_clones=20,dna_per_clone=200,A=abs(np.random.normal(loc=3,scale=1))+0.1+2,s=0.1,sigma_mu_ratio=30.,K=1.)
        RNA_sampled = np.array(RNA_sampled)/sum_RNA #Save as numpy array and count-normalize (make sure total counts (K) are held constant!)
        DNA_sampled = np.array(DNA_sampled)/sum_DNA
        
        RNA_data2 = np.append(RNA_data2,RNA_sampled)
        DNA_data2 = np.append(DNA_data2,DNA_sampled)
        sequence_indicators2 = np.append(sequence_indicators2,sequence_number*np.ones(len(RNA_sampled)))
    all_zeros = (RNA_data2==0)&(DNA_data2==0)
    RNA_data2 = RNA_data2[~all_zeros]
    DNA_data2 = DNA_data2[~all_zeros]
    sequence_indicators2 = sequence_indicators2[~all_zeros]

    
    data_DF = pd.DataFrame(data={"RNA_data1":RNA_data1,"DNA_data1":DNA_data1,"RNA_data2":RNA_data2,"DNA_data2":DNA_data2,"sequence_indicators1":sequence_indicators1,"sequence_indicators2":sequence_indicators2})
    b2_mean, b2_std = MPRAudit_Functions.MPRAudit_function(data_DF,ratiofunction=1,paired=True,timepoints=1,numtrials=100,jackpow=3./5)
    print("Total time = "+repr(process_time()-starttime))
    print("b2_mean = "+repr(b2_mean))
    print("b2_std = "+repr(b2_std))
        

