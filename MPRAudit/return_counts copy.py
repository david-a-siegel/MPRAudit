import numpy as np
import scipy.stats as st

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