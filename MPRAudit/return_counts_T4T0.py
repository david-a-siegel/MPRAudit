import numpy as np
import scipy.stats as st

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

