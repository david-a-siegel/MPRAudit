import numpy as np
import os
import scipy.stats as st
import pandas as pd
# import argparse
from MPRAudit_functions import *

def main():

    #Generate some random data with variable "A", the RNA/DNA ratio, to create nonzero b^2:    
    RNA_data = np.array([])
    DNA_data = np.array([])
    sequence_indicators = np.array([])
    for sequence_number in range(5):
        RNA_sampled, DNA_sampled, sum_RNA, sum_DNA = return_counts(num_clones=20,dna_per_clone=200,A=abs(np.random.normal(loc=3,scale=2))+0.1,s=0.1,sigma_mu_ratio=30.,K=1.)
        RNA_sampled = np.array(RNA_sampled)/sum_RNA #Save as numpy array and count-normalize (make sure total counts (K) are held constant!)
        DNA_sampled = np.array(DNA_sampled)/sum_DNA
        
        RNA_data = np.append(RNA_data,RNA_sampled)
        DNA_data = np.append(DNA_data,DNA_sampled)
        sequence_indicators = np.append(sequence_indicators,sequence_number*np.ones(len(RNA_sampled)))
    
    data_DF = pd.DataFrame(data={"RNA_data":RNA_data,"DNA_data":DNA_data,"sequence_indicators":sequence_indicators})
    b2_mean, b2_std = MPRAudit_function(data_DF,ratiofunction=1,paired=False,timepoints=1,numtrials=10,jackpow=3./5)
    print("b2_mean = "+repr(b2_mean))
    print("b2_std = "+repr(b2_std))
                            

if __name__ == '__main__':
     main()

