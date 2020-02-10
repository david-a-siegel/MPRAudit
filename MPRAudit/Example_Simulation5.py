import numpy as np
import os
import scipy.stats as st
import pandas as pd
import MPRAudit_Functions
from time import process_time

def main():

    starttime = process_time()
    #Generate some random data with variable "A", the RNA/DNA ratio, to create nonzero b^2:    
    RNA_data1_T0 = np.array([])
    DNA_data1_T0 = np.array([])
    RNA_data1_T4 = np.array([])
    DNA_data1_T4 = np.array([])
    RNA_data2_T0 = np.array([])
    DNA_data2_T0 = np.array([])
    RNA_data2_T4 = np.array([])
    DNA_data2_T4 = np.array([])

    sequence_indicators1_T0 = np.array([])
    sequence_indicators1_T4 = np.array([])
    sequence_indicators2_T0 = np.array([])
    sequence_indicators2_T4 = np.array([])

    #First do dataset 1, then do dataset 2 and change the T4T0_ratio
    for sequence_number in range(500):
        RNA_sampled_T0, DNA_sampled_T0, RNA_sampled_T4, DNA_sampled_T4, sum_RNA_T0, sum_DNA_T0, sum_RNA_T4, sum_DNA_T4 = MPRAudit_Functions.return_counts_T4T0(num_clones=20,dna_per_clone=200,A=abs(np.random.normal(loc=3,scale=.1))+0.1,s=0.1,sigma_mu_ratio=30.,T4T0_ratio=5,K=10.)
        RNA_sampled1_T0 = np.array(RNA_sampled_T0)/sum_RNA_T0 #Save as numpy array and count-normalize (make sure total counts (K) are held constant!)
        DNA_sampled1_T0 = np.array(DNA_sampled_T0)/sum_DNA_T0
        RNA_sampled1_T4 = np.array(RNA_sampled_T4)/sum_RNA_T4
        DNA_sampled1_T4 = np.array(DNA_sampled_T4)/sum_DNA_T4
                
        RNA_data1_T0 = np.append(RNA_data1_T0,RNA_sampled1_T0)
        DNA_data1_T0 = np.append(DNA_data1_T0,DNA_sampled1_T0)
        RNA_data1_T4 = np.append(RNA_data1_T4,RNA_sampled1_T4)
        DNA_data1_T4 = np.append(DNA_data1_T4,DNA_sampled1_T4)


        RNA_sampled_T0, DNA_sampled_T0, RNA_sampled_T4, DNA_sampled_T4, sum_RNA_T0, sum_DNA_T0, sum_RNA_T4, sum_DNA_T4 = MPRAudit_Functions.return_counts_T4T0(num_clones=20,dna_per_clone=200,A=abs(np.random.normal(loc=3,scale=.1))+0.1,s=0.1,sigma_mu_ratio=30.,T4T0_ratio=1,K=10.)
        RNA_sampled2_T0 = np.array(RNA_sampled_T0)/sum_RNA_T0 #Save as numpy array and count-normalize (make sure total counts (K) are held constant!)
        DNA_sampled2_T0 = np.array(DNA_sampled_T0)/sum_DNA_T0
        RNA_sampled2_T4 = np.array(RNA_sampled_T4)/sum_RNA_T4
        DNA_sampled2_T4 = np.array(DNA_sampled_T4)/sum_DNA_T4
                
        RNA_data2_T0 = np.append(RNA_data2_T0,RNA_sampled2_T0)
        DNA_data2_T0 = np.append(DNA_data2_T0,DNA_sampled2_T0)
        RNA_data2_T4 = np.append(RNA_data2_T4,RNA_sampled2_T4)
        DNA_data2_T4 = np.append(DNA_data2_T4,DNA_sampled2_T4)



        sequence_indicators1_T0 = np.append(sequence_indicators1_T0,sequence_number*np.ones(len(RNA_sampled_T0))) #This is a shortcut, not generally going to work perfectly
        sequence_indicators1_T4 = np.append(sequence_indicators1_T4,sequence_number*np.ones(len(RNA_sampled_T4)))
        sequence_indicators2_T0 = np.append(sequence_indicators2_T0,sequence_number*np.ones(len(RNA_sampled_T0)))
        sequence_indicators2_T4 = np.append(sequence_indicators2_T4,sequence_number*np.ones(len(RNA_sampled_T4)))

#     all_zeros_T0 = (RNA_data_T0==0)&(DNA_data_T0==0)  ## Let's just simulate data without too many zeros for now -- in general might want to remove clones with no counts at all
#     all_zeros_T4 = (RNA_data_T4==0)&(DNA_data_T4==0)
#     RNA_data_T0 = RNA_data_T0[~all_zeros_T0]
#     DNA_data_T0 = DNA_data_T0[~all_zeros_T0]
#     RNA_data_T4 = RNA_data_T4[~all_zeros_T4]
#     DNA_data_T4 = DNA_data_T4[~all_zeros_T4]
# 
#     sequence_indicators_T0 = sequence_indicators_T0[~all_zeros_T0]
#     sequence_indicators_T4 = sequence_indicators_T4[~all_zeros_T4]

    
    data_DF = pd.DataFrame(data={"RNA_data1_T0":RNA_data1_T0,"DNA_data1_T0":DNA_data1_T0,"RNA_data1_T4":RNA_data1_T4,"DNA_data1_T4":DNA_data1_T4,"RNA_data2_T0":RNA_data2_T0,"DNA_data2_T0":DNA_data2_T0,"RNA_data2_T4":RNA_data2_T4,"DNA_data2_T4":DNA_data2_T4,"sequence_indicators1_T0":sequence_indicators1_T0,"sequence_indicators1_T4":sequence_indicators1_T4,"sequence_indicators2_T0":sequence_indicators2_T0,"sequence_indicators2_T4":sequence_indicators2_T4})
    b2_mean, b2_std = MPRAudit_Functions.MPRAudit_function(data_DF,ratiofunction=11,paired=True,timepoints=2,numtrials=100,jackpow=3./5)
    print("Total time = "+repr(process_time()-starttime))
    print("b2_mean = "+repr(b2_mean))
    print("b2_std = "+repr(b2_std))
                            

if __name__ == '__main__':
     main()

