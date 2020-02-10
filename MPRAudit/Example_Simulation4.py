import numpy as np
import os
import scipy.stats as st
import pandas as pd
import MPRAudit_Functions
from time import process_time

def main():

    starttime = process_time()
    #Generate some random data with variable "A", the RNA/DNA ratio, to create nonzero b^2,
    #but only keep RNA counts (maybe it's a CRISPR screen):    
    RNA_data = np.array([])
    sequence_indicators = np.array([])

    for sequence_number in range(500):
        RNA_sampled, DNA_sampled, sum_RNA, sum_DNA = MPRAudit_Functions.return_counts(num_clones=20,dna_per_clone=100,A=abs(np.random.normal(loc=3,scale=.1))+0.1,s=0.1,sigma_mu_ratio=30.,K=10.)
        RNA_sampled = np.array(RNA_sampled) #Maybe we won't count-normalize this time               
        RNA_data = np.append(RNA_data,RNA_sampled)
        sequence_indicators = np.append(sequence_indicators,sequence_number*np.ones(len(RNA_sampled)))
    
    data_DF = pd.DataFrame(data={"RNA_data":RNA_data,"sequence_indicators":sequence_indicators})
    b2_mean, b2_std = MPRAudit_Functions.MPRAudit_function(data_DF,CRISPR_log_flag=True,paired=False,numtrials=100,jackpow=3./5)
    print("Total time = "+repr(process_time()-starttime))
    print("b2_mean = "+repr(b2_mean))
    print("b2_std = "+repr(b2_std))
                            

if __name__ == '__main__':
     main()

